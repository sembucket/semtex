///////////////////////////////////////////////////////////////////////////////
// normalise.C: scale a flow field so that its 2-norm is 1. Optionally
// select scale factor. This utility can be used to normalise
// eigenvectors.
//
// Copyright (c) 2002 <--> $Date$, Hugh Blackburn
//
// USAGE
// -----
// normalise [-h] [-s <num>] file
//
// NB: The 2-norm is not the same as the L2-norm, as the 2-norm does
// not account for the mass-matrix weighting of different nodal
// values.
//
// NBB: The input field is assumed to contain velocities and pressure ONLY!
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <cstring>

#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>

using namespace std;

#include <cfemdef.h>
#include <utility.h>
#include <blas.h>
#include <lapack.h>
#include <veclib.h>
#include <femlib.h>

static char prog[] = "normalise";

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

static void getargs   (int, char**, real_t&, ifstream&);
static void gethead   (istream&, hdr_info&);
static void allocate  (hdr_info&, vector<real_t*>&);
static void readdata  (hdr_info&, istream&, vector<real_t*>&);
static void scaledata (hdr_info&, const real_t&, vector<real_t*>&);
static void writedata (hdr_info&, ostream&, vector<real_t*>&);
static bool doswap    (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  ifstream        file;
  hdr_info        head;
  real_t          factor = 0.0;
  vector<real_t*> u;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, factor, file);
  gethead (file, head);

  allocate (head, u);
  readdata (head, file, u);

  file.close();
  
  scaledata (head, factor, u);
  writedata (head, cout, u);
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc  ,
		     char**    argv  ,
		     real_t&   factor,
		     ifstream& file  )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: normalise [-h] [-s <num>] file\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 's':
      if (*++argv[0]) factor  = atof (  *argv);
      else { --argc;  factor  = atof (*++argv); }
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1)
    file.open (argv[0], ios::in);
  else {
    cerr << usage;
    exit (EXIT_FAILURE);
  }

  if (!file) {
    cerr << prog << ": unable to open file" << endl;
    exit (EXIT_FAILURE);
  }
}


static void gethead (istream&  file,
		     hdr_info& head)
// ---------------------------------------------------------------------------
// Load data structure from file header info. Note that the field
// names are packed into a string without spaces.
// ---------------------------------------------------------------------------
{
  char  buf[StrMax];
  int_t i, j; 

  file.get (head.session, 25); file.getline (buf, StrMax);
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);

  file.get (head.created, 25); file.ignore (StrMax, '\n');

  file >> head.nr >> head.ns >> head.nz >> head.nel;
  file.ignore (StrMax, '\n');

  file >> head.step;     file.ignore (StrMax, '\n');
  file >> head.time;     file.ignore (StrMax, '\n');
  file >> head.timestep; file.ignore (StrMax, '\n');
  file >> head.kinvis;   file.ignore (StrMax, '\n');
  file >> head.beta;     file.ignore (StrMax, '\n');

  file.get (buf, 25);
  for (i = 0, j = 0; i < 25; i++)
    if (isalpha (buf[i])) head.fields[j++] = buf[i];
  file.ignore (StrMax, '\n');
  head.fields[j] = '\0';

  file.get (head.format, 25); file.getline (buf, StrMax);

  if (!strstr (head.format, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (head.format, "-endia"))
    message (prog, "input field file in unknown binary format", WARNING);
}


static void allocate (hdr_info&        head,
		      vector<real_t*>& u   )
// ---------------------------------------------------------------------------
// Allocate enough storage to hold all the data fields (sem format).
// ---------------------------------------------------------------------------
{
  const int_t nfield = strlen (head.fields);
  const int_t ntot   = head.nr * head.ns * head.nel * head.nz;
  int_t       i;

  u.resize (nfield);
  
  for (i = 0; i < nfield; i++) {
    u[i] = new real_t [ntot];
    Veclib::zero (ntot, u[i], 1);
  }
}


static void readdata (hdr_info&        head,
		      istream&         file,
		      vector<real_t*>& u   )
// ---------------------------------------------------------------------------
// Binary read of data areas, with byte-swapping if required.
// ---------------------------------------------------------------------------
{
  int_t   i, len;
  real_t* addr;
  bool    swab;

  const int_t nfield = strlen (head.fields);
  const int_t ntot   = head.nr * head.ns * head.nel * head.nz;

  // -- Read the base flow into the first plane location of u.

  swab = doswap (head.format);
  
  for (i = 0; i < nfield; i++) {
    addr = u[i];
    len  = ntot * sizeof (real_t);
    file.read (reinterpret_cast<char*>(addr), len);
    if (swab) Veclib::brev (ntot, addr, 1, addr, 1);
  }
}


static void scaledata (hdr_info&        head  ,
		       const real_t&    factor,
		       vector<real_t*>& u    )
// ---------------------------------------------------------------------------
// Scale the data such that the 2-norm of the velocities is unity
// (optionally a forced rescaling).
// ---------------------------------------------------------------------------
{
  int_t  i;
  real_t norm = 0.0;

  const int_t nfield = strlen (head.fields);
  const int_t ncomps = nfield - 1;
  const int_t ntot   = head.nr * head.ns * head.nel * head.nz;

  if (factor == 0.0) {
    for (i = 0; i < ncomps; i++) norm += Blas::nrm2 (ntot, u[i], 1);
    for (i = 0; i < nfield; i++) Blas::scal (ntot, 1.0/norm, u[i], 1);
  } else
    for (i = 0; i < nfield; i++) Blas::scal (ntot, factor,   u[i], 1);
}


static void writedata (hdr_info&        head,
		       ostream&         file,
		       vector<real_t*>& u   )
// ---------------------------------------------------------------------------
// Write out the data, semtex format.
// ---------------------------------------------------------------------------
{
  const int_t nfield = strlen (head.fields);
  const int_t ntot   = head.nr * head.ns * head.nel * head.nz;
  char        buf[StrMax], tmp[StrMax];
  int_t       i, j;

  sprintf (buf, hdr_fmt[0], head.session);
  file << buf;
  sprintf (buf, hdr_fmt[1], head.created);
  file << buf;

  sprintf (tmp, "%1d %1d %1d %1d", head.nr, head.ns, head.nz, head.nel);
  sprintf (buf, hdr_fmt[2], tmp);
  file << buf;

  sprintf (buf, hdr_fmt[3], head.step);     file << buf;
  sprintf (buf, hdr_fmt[4], head.time);     file << buf;
  sprintf (buf, hdr_fmt[5], head.timestep); file << buf;
  sprintf (buf, hdr_fmt[6], head.kinvis);   file << buf;
  sprintf (buf, hdr_fmt[7], head.beta);     file << buf;

  sprintf (buf, hdr_fmt[8], head.fields);
  file << buf;  

  sprintf (tmp, "binary ");
  Veclib::describeFormat (tmp + strlen (tmp));
  sprintf (buf, hdr_fmt[9], tmp);
  file << buf;

  for (i = 0; i < nfield; i++) 
    file.write (reinterpret_cast<char*>(u[i]), ntot * sizeof (real_t));
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
