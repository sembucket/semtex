///////////////////////////////////////////////////////////////////////////////
// flipmap.C: generate a list of index pairs for a symmetry-related transform.
//
// Copyright (C) Hugh Blackburn 2002.
//
// USAGE
// flipmap [options] [file]
// options:
//   -h       ... print this message
//   -x || -y ... reflection symmetry to be used to generate map
//
// FILES
// Input file is a semtex/prism ASCII mesh file, generated e.g. by meshpr.
// Output is ASCII, with three header lines followed by a list of index pairs:
//   8 8 1 108 NR NS NZ NEL  # -- matches input file header
//   x                       # -- reflection symmetry generator
//   202                     # -- NFLIP, no. of pairs to follow < NR*NS*NEL/2
//   11   4001
//   12   4000
//   13   3999
//   ... etc. (NFLIP lines)
//
// The mapping between the mesh points may not be unique; we take the
// first mapping found for each case (reflect positive->negative,
// *and* negative->positive). This way a single gather will do the
// reflection, and leave no holes.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

static char prog[] = "flipmap";

const double EPS = EPSSP;	// Positional tolerance.

static void getargs  (int, char**, char&, ifstream&);
static int  header   (ifstream&);
static int  loadmesh (ifstream&, const int, const char,
		      vector<real>&, vector<real>&);
static void findmap  (const char,
		      const int, const vector<real>&, const vector<real>&,
		      const int, vector<int>&, vector<int>& );
static void printup  (const char, const int,
		      const vector<int>&, const vector<int>&);


int main (int    argc,
          char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  ifstream     file;
  char         generator;
  vector<real> x, y;
  vector<int>  orig, flip;
  int          npts, nmap;

  getargs (argc, argv, generator, file);

  npts = header   (file);
  nmap = loadmesh (file, npts, generator, x, y);
  
  findmap (generator, npts, x, y, nmap, orig, flip);
  printup (generator, nmap, orig, flip);

  return EXIT_SUCCESS;
}


static void getargs (int       argc,
		     char**    argv,
		     char&     gen ,
		     ifstream& file)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: flipmap [-h] -x || -y [meshfile]\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'x':
      gen = 'x';
      break;
    case 'y':
      gen = 'y';
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (!(gen == 'x' || gen == 'y')) {
    cerr << prog << ": must specify -x or -y" << endl;
    exit (EXIT_FAILURE);
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


static int header (ifstream& file)
// ---------------------------------------------------------------------------
// If file header indicates this is a mesh file, output header and
// return number of points. Else die.
// ---------------------------------------------------------------------------
{
  int  nr, ns, nz, nel;
  char buf[StrMax];
  
  file >> nr >> ns >> nz >> nel;
  file.getline (buf, StrMax);

  if (!strstr (buf, "NR NS NZ NEL")) {
    cerr << prog
	 << ": mesh header line should include NR NS NZ NEL: " 
	 << buf << endl;
    exit (EXIT_FAILURE);
  }
  
  cout << nr << " " << ns << " " << 1 << " " << nel << " NR NS NZ NEL" << endl;

  return nr * ns * nel;
}


static int loadmesh (ifstream&     file,
		     const int     npts,
		     const char    gen ,
		     vector<real>& x   ,
		     vector<real>& y   )
// ---------------------------------------------------------------------------
// Load the vectors x and y. Return the number of points to be mapped
// by symmetry.
// ---------------------------------------------------------------------------
{
  int i, nmap = 0;

  x.setSize (npts);
  y.setSize (npts);

  for (i = 0; file && i < npts; i++) file >> x[i] >> y[i];

  if (!file) {
    cerr << prog << ": premature end of input" << endl;
    exit (EXIT_FAILURE);
  }

  if (gen == 'x') {
    for (i = 0; i < npts; i++) if (x[i] > EPS) nmap++;
  } else {
    for (i = 0; i < npts; i++) if (y[i] > EPS) nmap++;
  }

  nmap *= 2;

  if (nmap == 0) {
    cerr << prog << ": no points to map" << endl;
    exit (EXIT_FAILURE);
  }

  if (nmap > npts) {
    cerr << prog << ": too many points to map: "
	 << nmap << " vs. " << npts << endl;
    exit (EXIT_FAILURE);
  }

  return nmap;
}


static void findmap  (const char          gen ,
		      const int           npts,
		      const vector<real>& x   ,
		      const vector<real>& y   ,
		      const int           nmap,
		      vector<int>&        orig,
		      vector<int>&        flip)
// ---------------------------------------------------------------------------
// This is where the mapping gets constructed. Order npts*npts operation.
// Take the first available mapping index for each point.
// ---------------------------------------------------------------------------
{
  int i, j, k = 0, found;

  orig.setSize (nmap);
  flip.setSize (nmap);

  if (gen == 'x') {
    for (i = 0; i < npts; i++) {
      if (x[i] > EPS || x[i] < -EPS) {
	for (found = 0, j = 0; !found && j < npts; j++) {
	  if (fabs(x[i] + x[j]) < EPS && fabs(y[i] - y[j]) < EPS) {
	    orig[k]   = i;
	    flip[k++] = j;
	    found     = 1;
	  }
	}
	if (!found) {
	  cerr << "Warning: mesh point "
	       << x[i] << ",\t" << y[i] << "\tnot mirrored" << endl;
	}
      }
    }
  } else {
    for (i = 0; i < npts; i++) {
      if (y[i] > EPS || y[i] < -EPS) {
	for (found = 0, j = 0; !found && j < npts; j++) {
	  if (fabs(y[i] + y[j]) < EPS && fabs(x[i] - x[j]) < EPS) {
	    orig[k]   = i;
	    flip[k++] = j;
	    found     = 1;
	  }
	}
	if (!found) {
	  cerr << "Warning: mesh point "
	       << x[i] << ",\t" << y[i] << "\tnot mirrored" << endl;
	}
      } 
    }
  }

  if (k != nmap) {
    cerr << prog
	 << ": number of maps found, " << k
	 << ", not equal to number needed, " << nmap << endl;
    exit (EXIT_FAILURE);
  }
}


static void printup (const char         gen ,
		     const int          nmap,
		     const vector<int>& orig,
		     const vector<int>& flip)
// ---------------------------------------------------------------------------
// Print the map on standard output.
// ---------------------------------------------------------------------------
{
  int i;

  cout << gen  << endl;
  cout << nmap << endl;
  
  for (i = 0; i < nmap; i++)
    cout << orig[i] << '\t' << flip[i] << endl;
}

