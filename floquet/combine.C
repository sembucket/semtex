///////////////////////////////////////////////////////////////////////////////
// combine.C: add a perturbation mode to a base flow.
//
// Copyright (c) 2002 <--> $Date$, Hugh Blackburn
//
// USAGE
// -----
// combine [options] base pert
// options:
// -h       ... print this message.
// -b <num> ... set beta, wavenumber of base flow to <num> (3D) [Default: 1.0]
// -r <num> ... relative size of perturbation is <num>.         [Default: 1.0]
// -m <num> ... mode number for perturbation is <num> (3D only) [Default: 1]
// 
// Write to standard output.
//
// NOTES
// -----
// 1. Four cases of combinations of base and perturbation fields may
// be identified based on the number of velocity field variables and
// planes of data in the perturbation (see README file):
//
// N_BASE  N_PERT  N_Z    COMPLEXITY
//  2       2       1     Real only:    u.Re v.Re      p.Re   (O2_2D)
//  2       3       1     Half-complex: u.Re v.Re w.Im p.Re   (O2_3D_SYMM)
//  2       3       2     Full-complex (not implemented yet). (O2_3D)
//  3       3       2     Full-complex                        (SO2_3D)
//
// The structure of the output matches the above: whenever the
// perturbation is real only, so is the combined field (it has N_Z=1
// if N_BASE=2 and N_Z=2 if N_BASE=3); for all other cases the
// combined field has N_Z>=4 (i.e. is fully 3D).  For the half- and
// full-complex cases, the perturbation field is considered to be a
// Fourier mode.  The output field is transformed to physical space.
//
// 2. Note that only the first data plane of the base flow is read in,
// so N_Z=2 is effectively the same as N_Z=1 for the base flow: if all
// we have is a N_Z=1 base flow it can safely be projected to N_Z=2
// for use by this program. In future this should be fixed, but for
// now the workaround is OK.
//
// 3. Data files are assumed to only have velocities and pressure
// data, in the standard ordering, i.e. uvp or uvwp. The pressure of
// the outcome is set to zero, and is ignored in doing the scaling.
//
// 4. When producing a 3D field file, the file is given minimal number
// of Fourier modes required to contain it, subject to the 2, 3, 5
// prime factor requirement for FFT: the number of modes is rounded up
// to suit this.
//
// 5. NB: The relative size of the perturbation is based on the
// 2-norms of the respective fields, not the same as the L2-norms,
// since the 2-norm does not account for the mass-matrix weighting of
// different nodal values.  This means that a perturbation of constant
// energy will appear larger where the mesh density is lower. This
// should really be fixed but it means that we will need the session
// file as well.
//
// 6. With -r 0 the perturbation field is not scaled - it is output
// with the same energy it had on input.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <cstring>

#include <iostream>
#include <fstream>

#include <iomanip>
#include <vector>

using namespace std;

#include <cfemdef.h>
#include <utility.h>
#include <blas.h>
#include <lapack.h>
#include <veclib.h>
#include <femlib.h>

static char prog[] = "combine";

static char* hdr_fmt[] = { 
  "%-25s "    "Session\n",
  "%-25s "    "Created\n",
  "%-25s "    "Nr, Ns, Nz, Elements\n",
  "%-25d "    "Step\n",
  "%-25.6g "  "Time\n",
  "%-25.6g "  "Time step\n",
  "%-25.6g "  "Kinvis\n",
  "%-25.6g "  "Beta\n",
  "%-25s "    "Fields written\n",
  "%-25s "    "Format\n"
};

typedef struct hdr_data {
  char   session[StrMax];
  char   created[StrMax];
  int_t  nr, ns, nz, nel;
  int_t  step;
  real_t time;
  real_t timestep;
  real_t kinvis;
  real_t beta;
  char   fields[StrMax];
  char   format[StrMax];
} hdr_info;

static void  getargs   (int, char**, int_t&, real_t&, real_t&,
			ifstream&, ifstream&);
static void  gethead   (istream&, hdr_info&);
static bool  conform   (const hdr_info&, const hdr_info&);
static int_t roundup   (int_t&, int_t&, const hdr_info&, const hdr_info&);
static void  allocate  (hdr_info&, const int_t, vector<real_t*>&);
static void  readdata  (hdr_info&, istream&, hdr_info&, istream&,
			vector<real_t*>&, const real_t);
static void  packdata  (hdr_info&, const int_t, const int_t, vector<real_t*>&);
static void  transform (hdr_info&, const int_t, vector<real_t*>&, const int_t);
static void  writedata (hdr_info&, ostream&, const int_t,
			const real_t, vector<real_t*>&);
static bool  doswap    (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  ifstream        bFile, pFile;	// -- b ==> base, p ==> perturbation.
  hdr_info        bHead, pHead;
  int_t           nBase, nPert;
  int_t           mode = 1, nz;
  real_t          wght = 1.0, beta = 1.0;
  vector<real_t*> u;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, mode, wght, beta, bFile, pFile);
  gethead (bFile, bHead);
  gethead (pFile, pHead);

  conform (bHead, pHead);
  roundup (mode, nz, bHead, pHead);

  allocate (pHead, nz, u);
  readdata (bHead, bFile, pHead, pFile, u, wght);
  bFile.close(); pFile.close();

  if (nz > 1) {
    packdata  (pHead, mode, nz, u);
    transform (pHead, nz, u, INVERSE);
  }
  writedata (pHead, cout, nz, beta, u);
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     int_t&    mode ,
		     real_t&   wght ,
		     real_t&   beta ,
		     ifstream& bfile,
		     ifstream& pfile)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: combine [options] base pert\n"
    "options:\n"
    "-h       ... print this message.\n"
    "-b <num> ... set beta, wavenumber of base flow to <num> (3D)"
    " [Default: 1.0]\n"
    "-r <num> ... relative size of perturbation is <num>."
    " [Default: 1.0]\n"
    "-m <num> ... mode number for perturbation is <num> (3D only)"
    " [Default: 1]\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'b':
      if (*++argv[0])  beta = atof (*argv);
      else           { beta = atof (*++argv); argc--; }
      break;
    case 'm':
      if (*++argv[0])  mode = atoi (*argv);
      else           { mode = atoi (*++argv); argc--; }
      break;
    case 'r':
      if (*++argv[0])  wght = atof (*argv);
      else           { wght = atof (*++argv); argc--; }
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 2) {
    bfile.open (argv[0], ios::in);
    pfile.open (argv[1], ios::in);
  } else {
    cerr << usage;
    exit (EXIT_FAILURE);
  }

  if (!bfile) {
    cerr << prog << ": unable to open base file" << endl;
    exit (EXIT_FAILURE);
  }

  if (!pfile) {
    cerr << prog << ": unable to open perturbation file" << endl;
    exit (EXIT_FAILURE);
  }

  if (mode < 1) {
    cerr << prog << ": mode number must be a positive int" << endl;
    exit (EXIT_FAILURE);
  }

  if (wght < 0.0) {
    cerr << prog << ": perturbation weight must be positive" << endl;
    exit (EXIT_FAILURE);
  }
}


static void gethead (istream&  file  ,
		     hdr_info& header)
// ---------------------------------------------------------------------------
// Load data structure from file header info. Note that the field
// names are packed into a string without spaces.
// ---------------------------------------------------------------------------
{
  char  buf[StrMax];
  int_t i, j; 

  file.get (header.session, 25); file.getline (buf, StrMax);
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);

  file.get (header.created, 25); file.ignore (StrMax, '\n');

  file >> header.nr >> header.ns >> header.nz >> header.nel;
  file.ignore (StrMax, '\n');

  file >> header.step;     file.ignore (StrMax, '\n');
  file >> header.time;     file.ignore (StrMax, '\n');
  file >> header.timestep; file.ignore (StrMax, '\n');
  file >> header.kinvis;   file.ignore (StrMax, '\n');
  file >> header.beta;     file.ignore (StrMax, '\n');

  file.get (buf, 25);
  for (i = 0, j = 0; i < 25; i++)
    if (isalpha (buf[i])) header.fields[j++] = buf[i];
  file.ignore (StrMax, '\n');
  header.fields[j] = '\0';

  file.get (header.format, 25); file.getline (buf, StrMax);

  if (!strstr (header.format, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (header.format, "-endia"))
    message (prog, "input field file in unknown binary format", WARNING);
}


static bool conform (const hdr_info& bhead,
		     const hdr_info& phead)
// ---------------------------------------------------------------------------
// Check that the base and perturbation field conform for combination.
// This means they have the same element orders and number of
// elements, also that the field scalar component strings agree --
// except in case where base flow is 2D,nC and perturbation flow is
// 3D,(n+1)C, where 'w' is the only allowed additional scalar field,
// assumed to be the third.  Die if these requirements aren't
// satisfied.  Return a flag to indicate if base field requires
// padding for 'w'.
// ---------------------------------------------------------------------------
{
  bool pad = false;
  
  if (bhead.nr != phead.nr || bhead.ns != phead.ns || bhead.nel != phead.nel)
    message (prog, "base and perturbation sizes do not conform", ERROR);

  if (!((bhead.nz == 1 && phead.nz == 1) || (bhead.nz == 2 && phead.nz == 2)))
    message (prog, "3D structures of base and perturbation are bad", ERROR);

  if (strlen (bhead.fields) < 3 || strlen (phead.fields) < 3)
    message (prog, "base or perturbation don't have enough components", ERROR);

  if (bhead.nz == 1 && phead.nz == 1)
    if (!strcmp (bhead.fields, phead.fields))
      pad = false;
    else if (strlen (bhead.fields) == (strlen (phead.fields) - 1)) {
      int_t i, n = strlen (phead.fields);
      char extend[32];
      for (i = 0; i < 32; i++) extend[i] = '\0';
      extend[0] = bhead.fields[0];
      extend[1] = bhead.fields[1];
      extend[2] = 'w';
      for (i = 2; i < n; i++) extend[i+1] = bhead.fields[i];
      if (!strcmp (extend, phead.fields))
	pad = true;
      else
	message (prog, "can't extend base field list to perturbation", ERROR);
    }

  if (bhead.nz == 2 && phead.nz == 2)
    if (!strcmp (bhead.fields, phead.fields))
      pad = false;
    else
      message (prog, "base field list doesn't match perturbation", ERROR);

  return pad;
}


static int_t roundup (int_t&          mode ,
		      int_t&          nz   ,
		      const hdr_info& bhead,
		      const hdr_info& phead)
// ---------------------------------------------------------------------------
// Decide the number of z planes for the output field, based on
// restrictions outlined at the top of this file.  If the output field
// will be 3D, calculate nz from supplied mode number, then increment
// until it suits FFT.
// ---------------------------------------------------------------------------
{
  int_t n, ip, iq, ir, ipqr2;

  if (bhead.nz == 1        &&  phead.nz == 1        &&
      strlen(phead.fields) == strlen(bhead.fields)) {
    mode = 0;
    nz   = 1;
  } else {
    nz = 2 * (mode + 1);
    do {
      n = nz;
      Femlib::primes235 (n, ip, iq, ir, ipqr2);
      nz += 2;
    } while (n == 0);
    nz = n;
  }
  return nz;
}


static void allocate (hdr_info&        header,
		      const int_t      nz    ,
		      vector<real_t*>& u     )
// ---------------------------------------------------------------------------
// Allocate enough storage to hold all the data fields (sem format).
// Return the length of each scalar field (padded also so it's even).
// Initialise all the storage to zero.
// ---------------------------------------------------------------------------
{
  const int_t nfield    = strlen (header.fields);
  const int_t ntotelmt  = header.nr * header.ns * header.nel;
  const int_t planesize = ntotelmt + (ntotelmt % 2);
  const int_t ntot      = planesize * nz;
  int_t       i;

  u.resize (nfield);
  
  for (i = 0; i < nfield; i++) {
    u[i] = new real_t [ntot];
    Veclib::zero (ntot, u[i], 1);
  }
}


static void readdata (hdr_info&        bhead,
		      istream&         bfile,
		      hdr_info&        phead,
		      istream&         pfile,
		      vector<real_t*>& u    ,
		      const real_t     eps  )
// ---------------------------------------------------------------------------
// Binary read of data areas, with byte-swapping if required.  The
// initial storage of base and perturbation in u is the same as for
// "dual" files if the combination will eventually be 3D: the first
// plane's worth of data is the base flow, the second and third are
// the perturbation.  If the end product will be 2D, the perturbation
// is added to the base flow in the same (zeroth) plane.  In order to
// be able to scale the perturbation appropriately, we have to read it
// into temporary storage first.  We may also have to pad the base
// flow by adding an extra storage field ('w'); the storage area has
// enough space to cope.
//
// If input scale eps=0 then the perturbation is not scaled.
// ---------------------------------------------------------------------------
{
  bool    swab;
  int_t   i, j, len;
  real_t* addr;
  real_t  u2 = 0.0, U2 = 0.0;

  const int_t nBfield   = strlen (bhead.fields);
  const int_t nPfield   = strlen (phead.fields);

  const int_t nzB       = bhead.nz;
  const int_t nzP       = phead.nz;
  const int_t ntotelmt  = bhead.nr * bhead.ns * bhead.nel;
  const int_t planesize = ntotelmt + (ntotelmt % 2);

  vector<vector<real_t> > utmp (nPfield);

  // -- Read the base flow into the first plane location of u.

  swab = doswap (bhead.format);
  
  for (i = 0; i < nBfield; i++) {
    addr = u[i];
    len  = ntotelmt * sizeof (real_t);
    bfile.read (reinterpret_cast<char*>(addr), len);
    if (swab) Veclib::brev (planesize, addr, 1, addr, 1);
    bfile.ignore ((nzB - 1) * ntotelmt * sizeof (real_t));
  }

  // -- If required, move the 3rd and higher fields up by 1 to allow for 'w'.
  //    After this, both the base and perturbation fields have same number
  //    of components, nPfield.

  if (nPfield == nBfield + 1) {
    for (i = nBfield; i > 2; i--)
      Veclib::copy (planesize, u[i-1], 1, u[i], 1);
    Veclib::zero (planesize, u[2], 1);
  }
  Veclib::zero (planesize, u[nPfield-1], 1); // -- Zero the pressure.

  // -- Read perturbation flow into temporary storage.
 
  swab = doswap (phead.format);
  
  for (i = 0; i < nPfield; i++) {
    utmp[i].resize (nzP*planesize);
    for (j = 0; j < nzP; j++) {
      len  = ntotelmt * sizeof (real_t);
      addr = &utmp[i][j*planesize];
      pfile.read (reinterpret_cast<char*>(addr), len);
      if (swab) Veclib::brev (planesize, addr, 1, addr, 1);
    }
  }

  // -- If eps!=0, compute the 2-norms of base and perturbation flows,
  // scale perturbation.  If the base flow's energy is zero (i.e. we
  // just want to look at the perturbation field, have used a zero
  // base field), reset it to be 1.0.  NB: ignore the pressure fields.

  if (fabs (eps) > EPSDP) {
    for (i = 0; i < nPfield-1; i++) {
      U2 += Blas::nrm2 (ntotelmt, u[i], 1);
      for (j = 0; j < nzP; j++) {
	u2 += Blas::nrm2 (ntotelmt, &utmp[i][j*planesize], 1);
      }
    }
    
    U2 = (U2 < EPSDP) ? 1.0 : U2;
  
    for (i = 0; i < nPfield; i++)
      for (j = 0; j < nzP; j++)
	Blas::scal (ntotelmt, eps*U2/u2, &utmp[i][j*planesize], 1);
  }

  // -- Add the scaled perturbation to the base flow.

  if (nPfield == nBfield && nzP == 1) {	// -- Everything has 1 plane of data.
    for (i = 0; i < nPfield-1; i++)
      Veclib::vadd (ntotelmt, &utmp[i][0], 1, u[i], 1, u[i], 1);
  } else if (nPfield == nBfield + 1 && nzP == 1) { // -- half-complex pert.
    for (i = 0; i < nPfield-1; i++)
      if (phead.fields[i] != 'w')
	Veclib::copy (ntotelmt, &utmp[i][0], 1, u[i]+planesize, 1);
      else
	Veclib::copy (ntotelmt, &utmp[i][0], 1, u[i]+2*planesize, 1);
  } else if (nPfield == nBfield && nzP == 2) { // -- full-complex pert.
    for (i = 0; i < nPfield-1; i++) {
      Veclib::copy (ntotelmt, &utmp[i][0],         1, u[i]  +planesize, 1);
      Veclib::copy (ntotelmt, &utmp[i][planesize], 1, u[i]+2*planesize, 1);
    }
  }
}


static void packdata (hdr_info&        header,
		      const int_t      mode  ,
		      const int_t      nz    ,
		      vector<real_t*>& u     )
// ---------------------------------------------------------------------------
// This is only called if nz > 1.  Data in u are considered to be in
// Fourier-transformed state, and stored, dual-like, in the first 3
// planes of u.  Move planes around to place them in required spots.
// Zero other data.
// ---------------------------------------------------------------------------
{
  const int_t    nfield    = strlen (header.fields);
  const int_t    ntotelmt  = header.nr * header.ns * header.nel;
  const int_t    planesize = ntotelmt + (ntotelmt % 2);
  const int_t    nblock    = 2 * planesize;
  vector<real_t> tmp (2 * planesize);
  int_t          i;

  for (i = 0; i < nfield; i++) {
    Veclib::copy (nblock, u[i] + planesize, 1, &tmp[0], 1);
    Veclib::zero ((nz - 1) * planesize, u[i] + planesize, 1);
    Veclib::copy (nblock, &tmp[0], 1, u[i] + mode * nblock, 1);
  }
}


static void transform (hdr_info&        header,
		       const int_t      nz    , 
		       vector<real_t*>& u     , 
		       const int_t      dir   )
// ---------------------------------------------------------------------------
// (Inverse, here) Fourier transformation in z.
// ---------------------------------------------------------------------------
{
  const int_t nfield    = strlen (header.fields);
  const int_t ntotelmt  = header.nr * header.ns * header.nel;
  const int_t planesize = ntotelmt + (ntotelmt % 2);
  int_t       i;

  for (i = 0; i < nfield; i++) Femlib::DFTr (u[i], nz, planesize, dir);
}


static void writedata (hdr_info&        header,
		       ostream&         file  ,
		       const int_t      nz    ,
		       const real_t     beta  ,
		       vector<real_t*>& u     )
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  const int_t nfield    = strlen (header.fields);
  const int_t ntotelmt  = header.nr * header.ns * header.nel;
  const int_t planesize = ntotelmt + (ntotelmt % 2);
  char        buf[StrMax], tmp[StrMax];
  int_t       i, j;

  sprintf (buf, hdr_fmt[0], header.session);
  file << buf;
  sprintf (buf, hdr_fmt[1], header.created);
  file << buf;

  sprintf (tmp, "%1d %1d %1d %1d", header.nr, header.ns, nz, header.nel);
  sprintf (buf, hdr_fmt[2], tmp);
  file << buf;

  sprintf (buf, hdr_fmt[3], header.step);     file << buf;
  sprintf (buf, hdr_fmt[4], header.time);     file << buf;
  sprintf (buf, hdr_fmt[5], header.timestep); file << buf;
  sprintf (buf, hdr_fmt[6], header.kinvis);   file << buf;
  sprintf (buf, hdr_fmt[7], beta);            file << buf;

  sprintf (buf, hdr_fmt[8], header.fields);
  file << buf;  

  sprintf (tmp, "binary ");
  Veclib::describeFormat (tmp + strlen (tmp));
  sprintf (buf, hdr_fmt[9], tmp);
  file << buf;

  for (i = 0; i < nfield; i++) 
    for (j = 0; j < nz; j++)
      file.write ((char*) (u[i] + j * planesize), ntotelmt * sizeof (real_t));
}


static bool doswap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping of input is required to make sense of input.
// ---------------------------------------------------------------------------
{
  char mfmt[StrMax];

  Veclib::describeFormat (mfmt);   

  if (!strstr (ffmt, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (ffmt, "-endia"))
    message (prog, "input field file in unknown binary format", WARNING);

  return (strstr (ffmt, "big") && strstr (mfmt, "little")) || 
         (strstr (mfmt, "big") && strstr (ffmt, "little"));
}
