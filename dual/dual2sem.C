///////////////////////////////////////////////////////////////////////////////
// dual2sem.C: convert a dual 3D file to a semtex/nekton 3D file.
//
// Copyright (c) 2000 Hugh Blackburn
//
// USAGE
// -----
// dual2sem [options] [file]
// options:
// -h       ... print this message.
// -f <num> ... the mode number of the 3D mode is <num>. [Default: 1]
// -i       ... invert (sem --> dual).
// 
// If file is not present, read from standard input.  Write to
// standard output.
//
// A dual field file describes a 3D data set, but carries only two
// Fourier modes, the zeroth (just real data) and one more mode
// (complex), with a nominated mode number.  So only three planes of
// data are held, and in addition the single complex mode is held in
// complex form (its first plane is real, its second is imaginary).
//
// Dual2sem projects a dual field file onto a nekton type field file,
// with the minimal number of modes required to contain it (subject to
// the 2, 3, 5 prime factor requirements for FFT).  Dual2sem cannot
// determine from the dual field file the number of the fundamental 3D
// mode (the header for the file is standard semtex format), so it
// needs a command-line argument to set it.  Likewise the mode number
// is needed in the inverse operation (in order to select the
// fundamental mode number from the semtex/nekton file).
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include <iostream.h>
#include <fstream.h>
#include <strstream.h>
#include <iomanip.h>

#include <femdef.h>
#include <Array.h>
#include <Utility.h>
#include <Blas.h>
#include <Lapack.h>
#include <Veclib.h>
#include <Femlib.h>

static char prog[] = "dual2sem";

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
  int    nr, ns, nz, nel;
  int    step;
  double time;
  double timestep;
  double kinvis;
  double beta;
  char   fields[StrMax];
  char   format[StrMax];
} hdr_info;

static void getargs   (int, char**, int&, int&, ifstream&);
static void roundup   (const int, int&);
static void getheader (istream&, hdr_info&);
static void allocate  (hdr_info&, const int, vector<real*>&);
static void readdata  (hdr_info&, istream&, const int, vector<real*>&);
static void packdata  (hdr_info&, const int, const int, vector<real*>&,
		       const int);
static void transform (hdr_info&, const int, vector<real*>&, const int);
static void writedata (hdr_info&, ostream&, const int, vector<real*>&);
static int  doswap    (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  ifstream      file;
  hdr_info      header;
  int           dir = FORWARD, mode = 1, nz;
  vector<real*> u;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, dir, mode, file);

  if (dir == FORWARD) {		// -- Project a dual field file.

    roundup   (mode, nz);
    getheader (file, header);
    allocate  (header, nz, u);
    readdata  (header, file, 3, u);
    packdata  (header, mode, nz, u, INVERSE);
    transform (header, nz, u, INVERSE);
    writedata (header, cout, nz, u);

  } else {			// -- Restrict a sem field file.

    getheader (file, header);
    allocate  (header, header.nz, u);
    readdata  (header, file, header.nz, u);
    transform (header, header.nz, u, FORWARD);
    packdata  (header, mode, header.nz, u, FORWARD);
    writedata (header, cout, 3, u);

  }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc,
		     char**    argv,
		     int&      dir ,
		     int&      mode,
		     ifstream& file)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: dual2sem [options] [file]\n"
    "options:\n"
    "-h       ... print this message\n"
    "-f <num> ... the mode number of the 3D mode is <num>. [Default: 1]\n"
    "-i       ... invert (sem --> dual).\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      dir = INVERSE;
      break;
    case 'f':
      if (*++argv[0]) mode = atoi (*argv);
      else {mode = atoi (*++argv); argc--;}
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if   (argc == 1) file.open   (*argv, ios::in);
  else             file.attach (0);

  if (!file) {
    cerr << prog << ": unable to open input file" << endl;
    exit (EXIT_FAILURE);
  }

  if (mode < 1) {
    cerr << prog << ": mode number must be a positive integer" << endl;
    exit (EXIT_FAILURE);
  }
}


static void roundup (const int  mode, 
		     int&       nz  )
// ---------------------------------------------------------------------------
// Calculate nz from mode number, then increment until it suits FFT.
// ---------------------------------------------------------------------------
{
  int n, ip, iq, ir, ipqr2;

  nz = 2 * (mode + 1);

  do {
    n = nz;
    Femlib::primes235 (n, ip, iq, ir, ipqr2);
    nz += 2;
  } while (n == 0);
  nz = n;
}


static void getheader (istream&  file  ,
		       hdr_info& header)
// ---------------------------------------------------------------------------
// Load data structure from file header info.
// ---------------------------------------------------------------------------
{
  char    buf[StrMax];
  integer i, j; 

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


static void allocate (hdr_info&       header,
		      const int       nz    ,
		      vector<real*>&  u     )
// ---------------------------------------------------------------------------
// Allocate enough storage to hold all the data fields (sem format).
// ---------------------------------------------------------------------------
{
  const int nfield    = strlen (header.fields);
  const int ntotelmt  = header.nr * header.ns * header.nel;
  const int planesize = ntotelmt + (ntotelmt % 2);
  const int ntot      = planesize * nz;
  int       i;

  u.setSize (nfield);
  
  for (i = 0; i < nfield; i++) {
    u[i] = new real [ntot];
    Veclib::zero (ntot, u[i], 1);
  }
}


static void readdata (hdr_info&      header,
		      istream&       file  ,
		      const int      nplane,
		      vector<real*>& u     )
// ---------------------------------------------------------------------------
// Binary read of data area, followed by byte-swapping if required.
// ---------------------------------------------------------------------------
{
  const int nfield    = strlen (header.fields);
  const int ntotelmt  = header.nr * header.ns * header.nel;
  const int planesize = ntotelmt + (ntotelmt % 2);
  const int ntot      = planesize * nplane;
  const int swab      = doswap (header.format);
  int       i, j;
  
  for (i = 0; i < nfield; i++) {
    for (j = 0; j < nplane; j++) {
      file.read ((char*) (u[i] + j * planesize), ntotelmt * sizeof (real));
    }
    if (swab) Veclib::brev (ntot, u[i], 1, u[i], 1);
  }
}


static void packdata (hdr_info&      header,
		      const int      mode  ,
		      const int      nz    ,
		      vector<real*>& u     ,
		      const int      dir   )
// ---------------------------------------------------------------------------
// Data in u are considered to be in Fourier-transformed state.  Move
// planes around to place them in required spots.  Zero other data.
// ---------------------------------------------------------------------------
{
  const int    nfield    = strlen (header.fields);
  const int    ntotelmt  = header.nr * header.ns * header.nel;
  const int    planesize = ntotelmt + (ntotelmt % 2);
  const int    nblock    = 2 * planesize;
  vector<real> tmp (2 * planesize);
  int          i;

  if (dir == FORWARD)		// -- Pack down to 3 data planes.
    for (i = 0; i < nfield; i++) {
      Veclib::copy (nblock, u[i] + mode * nblock, 1, tmp(), 1);
      Veclib::zero ((nz - 1) * planesize, u[i] + planesize, 1);
      Veclib::copy (nblock, tmp(), 1, u[i] + planesize, 1);
    }
  else				// -- Reposition nominated mode.
    for (i = 0; i < nfield; i++) {
      Veclib::copy (nblock, u[i] + planesize, 1, tmp(), 1);
      Veclib::zero ((nz - 1) * planesize, u[i] + planesize, 1);
      Veclib::copy (nblock, tmp(), 1, u[i] + mode * nblock, 1);
    }
}


static void transform (hdr_info&      header,
		       const int      nz    , 
		       vector<real*>& u     , 
		       const int      dir   )
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  const int nfield    = strlen (header.fields);
  const int ntotelmt  = header.nr * header.ns * header.nel;
  const int planesize = ntotelmt + (ntotelmt % 2);
  int       i;

  for (i = 0; i < nfield; i++) Femlib::DFTr (u[i], nz, planesize, dir);
}


static void writedata (hdr_info&      header,
		       ostream&       file  ,
		       const int      nz    ,
		       vector<real*>& u     )
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  const int nfield    = strlen (header.fields);
  const int ntotelmt  = header.nr * header.ns * header.nel;
  const int planesize = ntotelmt + (ntotelmt % 2);
  char      buf[StrMax], tmp[StrMax];
  int       i, j;

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
  sprintf (buf, hdr_fmt[7], header.beta);     file << buf;

  sprintf (buf, hdr_fmt[8], header.fields);
  file << buf;  

  sprintf (tmp, "binary ");
  Veclib::describeFormat (tmp + strlen (tmp));
  sprintf (buf, hdr_fmt[9], tmp);
  file << buf;

  for (i = 0; i < nfield; i++) 
    for (j = 0; j < nz; j++)
      file.write ((char*) (u[i] + j * planesize), ntotelmt * sizeof (real));
}


static int doswap (const char* ffmt)
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
