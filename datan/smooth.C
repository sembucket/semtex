///////////////////////////////////////////////////////////////////////////////
// smooth.cc: given a list of (equally-spaced) points, construct and apply
// a Savitsky--Golay FIR filter to generate smoothed output.  Optionally
// compute instead a derivative of desired order.  Symmetric (acausal) only.
//
// Usage: smooth [options] [file]
// options are:
// -d <num> ... compute derivative of desired order (default: 0)
// -s <num> ... stepsize, only used for derivatives (default: 1)
// -p <num> ... filter polynomial order             (default: 0)
// -w <num> ... filter semi-width                   (default: 1)
// -r       ... reproduce input in output
// -h       ... print this message
// -v       ... be verbose
//
// Note that the polynomial order must be less than (2 * width + 1).
// If not set, the width will be calculated to satisfy this criterion.
// Press et al suggest the order should be at least 4 for derivatives. 
//
// References: 
// Numerical Recipes, 2e, \S\,14.8.  Press et al, 1992.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

#include <Utility.h>
#include <Veclib.h>
#include <Array.h>
#include "nr77.h"

typedef double real;

static char prog[] = "smooth";
static void getargs (int,char**,int&,int&,int&,real&,int&,int&,istream*&);
static void filter  (istream&,ostream&,vector<real>&,const int&,const real&); 


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Drive routines in this file.
// ---------------------------------------------------------------------------
{
  istream*     file;
  char         err[StrMax];
  int          verbose = 0, repro = 0, npol = 0, nwid = 0, nder = 0;
  real         h = 1.0;
  vector<real> c;

  getargs (argc, argv, verbose, repro, nder, h, npol, nwid, file);

  if (nwid) {
    if (nwid < (npol + 1) >> 1) {
      sprintf (err, "semiwidth (%1d) must be at least %1d for order %1d",
	       nwid, (npol + 1) >> 1, npol);
      message (prog, err, ERROR);
    }
  } else
    nwid = (npol + 1) >> 1;

  c.setSize (2 * nwid + 1);
  Recipes::savgol (c(), c.getSize(), nwid, nwid, nder, npol);
  
  filter (*file, cout, c, nwid, 1.0 / pow (h, nder));

  return EXIT_SUCCESS;
}


static void getargs (int       argc   ,
		     char**    argv   ,
		     int&      verbose,
		     int&      repro  ,
		     int&      nder   ,
		     real&     h      ,
		     int&      npol   ,
		     int&      nwid   ,
		     istream*& file   )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char buf[StrMax], c;
  char usage[] =
    "Usage: %s [options] [file]\n"
    "options are:\n"
    "-d <num> ... compute derivative of desired order (default: 0)\n"
    "-s <num> ... stepsize, only used for derivatives (default: 1)\n"
    "-p <num> ... filter polynomial order             (default: 0)\n"
    "-w <num> ... filter semi-width                   (default: 1)\n"
    "-r       ... reproduce input in output\n"
    "-h       ... print this message\n"
    "-v       ... be verbose\n";
 
  while (--argc  && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'd':
      if   (*++argv[0])           nder = atoi (*argv);
      else              { --argc; nder = atoi (*++argv); }
      break;
    case 's':
      if   (*++argv[0])           h    = atof (*argv);
      else              { --argc; h    = atof (*++argv); }
      break;
    case 'p':
      if   (*++argv[0])           npol = atoi (*argv);
      else              { --argc; npol = atoi (*++argv); }
      break;
    case 'w':
      if   (*++argv[0])           nwid = atoi (*argv);
      else              { --argc; nwid = atoi (*++argv); }
      break;
    case 'v':
      do verbose++; while (*++argv[0] == 'v');
      break;
    case 'r':
      repro = 1;
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }
  
  if (argc == 1) {
    file = new ifstream (*argv);
    if (file -> bad()) message (prog, "unable to open input file", ERROR);
  } else file = &cin;
}


static void filter (istream&      ifile,
		    ostream&      ofile,
		    vector<real>& coeff,
		    const int&    N    ,
		    const real&   scale)
// ---------------------------------------------------------------------------
// Convolve input with filter, print on output.  Filter coefficients are
// supplied stored in coeff in wrap-around order.
// ---------------------------------------------------------------------------
{
  char         err[StrMax];
  register int i, j;
  real         z;
  const    int buflen = 2*N + 1;
  vector<real> buf (buflen);

  Veclib::zero (buflen, buf(), 1);
  
  // -- Read ahead to half-fill buf.

  for (i = 0; i < N; i++)
    ifile >> buf[N - i];
  if (!ifile) {
    sprintf (err, "couldn't find enough points to start up (needed %1d)", N);
    message (prog, err, ERROR);
  }

  // -- Deal with main body of input.

  while (ifile >> buf[0]) {
    z = buf[0] * coeff[0];
    for (i = 1; i <= N; i++) z += buf[i]*coeff[N + i] + buf[N + i]*coeff[i];
    for (i = buflen - 1; i > 0; i--) buf[i] = buf[i - 1];
    ofile << scale * z << endl;
  }

  // -- Clean out filter with zeros.

  for (j = 0; j < N; j++) {
    buf[0] = 0.0;
    z = buf[0] * coeff[0];
    for (i = 1; i <= N; i++) z += buf[i]*coeff[N + i] + buf[N + i]*coeff[i];
    for (i = buflen - 1; i > 0; i--) buf[i] = buf[i - 1];
    ofile << scale * z << endl;  
  }
}

 
