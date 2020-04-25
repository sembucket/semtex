///////////////////////////////////////////////////////////////////////////////
// meshplot.cpp: make a PostScript plot of a semtex/prism mesh from
// data in a file produced by the meshpr utility or similar.  Graphics
// are based on the PSplot package introduced in Numerical Recipes 3e
// Ch. 22 and described in webnote 26, http://wwww.nr.com/webnotes?26.
//
// Copyright (c) 2019 <--> $Date: 2019/05/30 06:36:12 $, Hugh Blackburn
//
// Usage:
// ------
// meshplot [options] [file]
// options:
//   -h        ... display this message
//   -a        ... do not show axes
//   -i        ... show element-internal structure
//   -n        ... show element numbers
//   -o 'file' ... write output to named file [Default: stdout]
//   -d 'prog' ... call prog to display named PostScript output file
//   -b 'xmin,xmax,ymin,ymax' ... limit view to region defined by string
//
// Files:
// ------
// Input consists of mesh information in ASCII format (produced by
// meshpr), with a single-line header followed by (x,y) mesh locations
// and optionally by z location list (which is ignored). File may be
// read on stdin.
//
// Mesh header is of form:
//
// 9 9 32 206 NR NS NZ NEL
//
// then follows (x,y) locations of each element's 2D mesh points in
// row-major order, ASCII format.  In all local implementations,
// NR=NS=np.  Finally the NZ + 1 z locations of the 3D mesh are
// supplied, if NZ>1 -- these are ignored here.  
//
// $Id: meshplot.cpp,v 9.1 2019/05/30 06:36:12 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>		// C standard headers.
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>

#include <vector>		// C++ standard headers.
#include <iostream>
#include <fstream>

#include <cfemdef.h>		// Semtex headers.
#include <utility.h>
#include <veclib.h>

using namespace std;

#include "psplot.h"

static char prog[] = "meshplot";
static void getargs  (int, char**, char*&, char*&, bool&, bool&, bool&,
		      double&, double&, double&, double&, istream*&);
static void readmesh (istream&, int&, int&, int&, int&,
		      vector<double>&, vector<double>&,
		      double&, double&, double&, double&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// An A4 page is 595 pts wide and 841 pts high, i.e. centred near 300, 420.
// ---------------------------------------------------------------------------
{
  char           *ofile = NULL, *disprog = NULL;
  FILE*          output;
  istream*       input;
  bool           internal = false, axes = true, numbers = false;
  int            i, j, k;
  int            nr, ns, nrns, nel, ntot;
  vector<double> x, y;
  double         xmin =  FLT_MAX, ymin =  FLT_MAX; // -- Obtained from mesh.
  double         xmax = -FLT_MAX, ymax = -FLT_MAX;
  double         Xmin =  FLT_MAX, Ymin =  FLT_MAX; // -- Using command-line.
  double         Xmax = -FLT_MAX, Ymax = -FLT_MAX;
  double         xlen, ylen, xcen, ycen;

  getargs (argc, argv, ofile, disprog, internal, axes, numbers,
	   Xmin, Xmax, Ymin, Ymax, input);

  readmesh (*input, nr, ns, nel, ntot, x, y, xmin, xmax, ymin, ymax);
  nrns = nr * ns;

  if (Xmin < Xmax) { // -- View limits were set on command line.
    xlen = Xmax - Xmin;
    ylen = Ymax - Ymin;
  } else {	    // -- Default: obtained from mesh file.
    xlen = xmax - xmin;
    ylen = ymax - ymin;
  }

  // -- Have all the info, now make the plot.

  double xlo, xhi, ylo, yhi;	// -- Limits on the page in pts.
                                // -- Maximum dimension 500 pts.

  if (xlen > ylen) {
    xlo = 300. - 250.;             xhi = 300. + 250.;
    ylo = 420. - ylen/xlen * 250.; yhi = 420. + ylen/xlen * 250.;
  } else {
    xlo = 300. - xlen/ylen * 250.; xhi = 300. + xlen/ylen * 250.;
    ylo = 420. - 250.;             yhi = 420. + 250.;
  }

  // -- Use psplot utilities (psplot.h) to do the plotting.

  output = (ofile) ? fopen (ofile, "w") : stdout;
  
  PSpage pg   (output);
  PSplot plot (pg, xlo, xhi, ylo, yhi);

  if (Xmin < Xmax) 
    plot.setlimits (Xmin, Xmax, Ymin, Ymax);
  else
    plot.setlimits (xmin, xmax, ymin, ymax);

  if (axes) {
    plot.frame();
    plot.autoscales();
  }

  plot.clip ();

  if (internal) {  // -- Element internal structure, draw first.
    plot.setgray      (0.5);
    plot.setlinewidth (0.5);

    for (i = 0; i < nel; i++) {
      for (j = 1; j < (ns - 1); j++) {
	// -- Lines of constant s.
	k = i * nrns + j * nr;
	plot.polyline (nr, &x[k],  1,  &y[k], 1);
      }
      for (j = 1; j < (nr - 1); j++) {
	// -- Lines of constant r.
	k = i * nrns + j;
	plot.polyline (nr, &x[k], -nr, &y[k], -nr);
      }
    }
  }

  // -- Element outlines.

  plot.setgray      (0.);
  plot.setlinewidth (1.);

  for (i = 0; i < nel; i++) {
    // -- Side 1.
    k = i * nrns;
    plot.polyline (nr, &x[k],  1,  &y[k], 1);
    // -- Side 2.
    k = i * nrns + (nr - 1);
    plot.polyline (ns, &x[k],  nr, &y[k], nr);
    // -- Side 3.
    k = i * nrns + (nr - 1) * ns;
    plot.polyline (nr, &x[k], -1,  &y[k], -1);
    // -- Side 4.
    k = i * nrns;
    plot.polyline (nr, &x[k], -nr, &y[k], -nr);
  }

  if (numbers) {  // -- Add element numbers at element centroids.
    char elmt[16];
    
    for (i = 0; i < nel; i++) {
      k = i * nrns;
      xcen = Veclib::sum (nrns, &x[k], 1) / nrns;
      ycen = Veclib::sum (nrns, &y[k], 1) / nrns;
      sprintf (elmt, "%1d", i+1);
      plot.Clabel (elmt, xcen, ycen);
    }
  }

  if (disprog && ofile) {
    fclose (output);
    plot.display (disprog, ofile);
  }

  return (EXIT_SUCCESS);
}


static void getargs (int       argc    ,
		     char**    argv    ,
		     char*&    ofile   ,
		     char*&    disprog ,
		     bool&     internal,
		     bool&     axes    ,
		     bool&     numbers ,
		     double&   Xmin    ,
		     double&   Xmax    ,
		     double&   Ymin    ,
		     double&   Ymax    ,
		     istream*& input   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: meshplot [options] [file]\n"
    "  options:\n"
    "  -h        ... display this message\n"
    "  -a        ... do not show axes\n"
    "  -i        ... show element-internal mesh\n"
    "  -n        ... show element numbers\n"
    "  -o <file> ... write output to named file [Default: stdout]\n"
    "  -d <prog> ... call prog to display named PostScript output file\n"
    "  -b 'xmin,xmax,ymin,ymax' ... limit view to region defined by string\n";
  char err[StrMax];
  char *tok, *pspec;
  int  set = 0;
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'a':
      axes = false;
      break;
    case 'i':
      internal = true;
      break;
    case 'n':
      numbers = true;
      break;
    case 'd':
      --argc;
      disprog = *++argv;
      break;
    case 'o':
      --argc;
      ofile = *++argv;
      break;
    case 'b':
      if (*++argv[0])
	pspec = *argv;
      else {
	--argc;
	pspec = *++argv;
      }
      if (tok = strtok (pspec, ",")) { 
	Xmin = atof (tok);
	set = 1;
      } else
	message (prog, "couldn't parse xmin from box string", ERROR);
      while (tok = strtok (0, ","))
	switch (++set) {
	case 2:
	  Xmax = atof (tok);
	  break;
	case 3:
	  Ymin = atof (tok);
	  break;
	case 4:
	  Ymax = atof (tok);
	  break;
	default:
	  message (prog, "too many numbers in box string", ERROR);
	    break;
	  }
	if (set != 4) {
	  sprintf (err, "wrong number of parameters in box string (%1d)",set);
	  message (prog, err, ERROR);
	}
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> fail()) message (prog, "unable to open geometry file", ERROR);
  } else input = &cin;
}


static void readmesh (istream&        file,
		      int&            nr  ,
		      int&            ns  ,
		      int&            nel ,
		      int&            ntot,
		      vector<double>& x   ,
		      vector<double>& y   ,
		      double&         xmin,
		      double&         xmax,
		      double&         ymin,
		      double&         ymax)
// ---------------------------------------------------------------------------
// File is already open; extract mesh info.
// ---------------------------------------------------------------------------
{
  const char routine[] = "readmesh";
  char  buf[STR_MAX], err[STR_MAX];
  int   i, nz;
 
  file >> nr >> ns >> nz >> nel;
  file.getline (buf, StrMax);

  if (!strstr (buf, "NR NS NZ NEL")) {
    sprintf (err, "mesh header line should include NR NS NZ NEL: %s", buf);
    message (routine, err, ERROR);
  }

  ntot = nr * ns * nel;
  x.resize (ntot);
  y.resize (ntot);

  for (i = 0; i < ntot; i++) {
    file >> x[i] >> y[i];
    xmin = min (xmin, x[i]);
    xmax = max (xmax, x[i]);
    ymin = min (ymin, y[i]);
    ymax = max (ymax, y[i]);
  }

  if (file.fail()) message (prog, "problem reading mesh", ERROR);
}

