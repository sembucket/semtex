///////////////////////////////////////////////////////////////////////////////
// cycint.C: compute definite integrals of cycles of a waveform.
// waveform.
//
// Usage:
// ------
// cycint [-h] [-v] [file]
//
// Synopsis:
// ---------
// Expected input is 3 columns of ASCII, the first is time, the second is
// a process whose zero crossings are used to define cycle times of the 
// process in the third column.  Find and print local and average values
// of the cyclic integrals of the third waveform, using cubic splines.
//
// Options:
// --------
// -h (help):    print usage string.
// -v (verbose): peaks, local and running average frequencies.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <strstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include <Utility.h>
#include <Stack.h>
#include <Array.h>
#include "nr77.h"

class triplet {
public:
  triplet (double X, double Y, double Z) : x (X), y (Y), z (Z) { }
  double    x, y, z;
};

static char prog[] = "cycint";
static int     num;
static double  sign;
static double* xx;
static double* yy;
static double* cc;
static double  wave     (const double&);
static int     dcycle   (const int, const double*, const double*, 
		         const double*, double&, double&);

extern "C" {
double         dsplquad (const double*, const double*, const double*,
			 const int, const double, const double);
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  triplet*        datum;
  Stack<triplet*> data;
  vector<double>  t, r, f;	// -- Time, reference, function.
  int             i, N;

  double t1, t2, r1, r2, f1, f2, period, end;
  double freq, favg = 0.0;
  double area, aavg = 0.0;

  char   line[StrMax], c;
  int    crossed = 0, verbose = 0, ncycles = 0;
  char   usage[] = "cycint [-h] [-v] [file]";

  // -- Process command line.

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      cout << usage << endl;
      return EXIT_SUCCESS;
      break;
    case 'v':
      verbose = 1;
      break;
    default:
      cerr << usage << endl;
      return EXIT_FAILURE;
      break;
    }

  if (argc == 1) {
    ifstream* inputfile = new ifstream (*argv);
    if (inputfile -> good()) {
      cin = *inputfile;
      } else {
	cerr << prog << "unable to open file" << endl;
	exit (EXIT_FAILURE);
      }
  }

  // -- Find first crossing.

  cin >> t1 >> r1 >> f1;
  while (!crossed && cin >> t2 >> r2 >> f2) {
    if (crossed = ((r1 < 0.0) && (r2 >= 0.0))) {
      datum = new triplet (t1, r1, f1);
      data.push (datum);
      if (verbose) {
	end = t1 + (t1 - t2) * r1 / (r2 - r1);
	sprintf (line, "%.4e", end);
	cout << line << endl;
      }
    } 
    t1 = t2;
    r1 = r2;
    f1 = f2;
  }

  if (!crossed) {
    cerr << "cycint: no crossing found" << endl;
    return EXIT_FAILURE;

  } else {
    datum = new triplet (t2, r2, f2);
    data.push (datum);
    crossed = 0;
  }

  // -- Do the rest of input, cycle by cycle.

  while (cin) {

    while (!crossed && cin >> t2 >> r2 >> f2) {
      crossed = ((r1 < 0.0) && (r2 >= 0.0));
      datum = new triplet (t2, r2, f2);
      data.push (datum);
      t1 = t2;
      r1 = r2;
      f1 = f2;
    }

    if (crossed) {		// -- End of current cycle.
      i = N = data.depth();
      t.setSize (N);
      r.setSize (N);
      f.setSize (N);

      ncycles += 1;
      crossed  = 0;

      while (i--) {
	datum = data.pop();
	t[i]  = datum -> x;
	r[i]  = datum -> y;
	f[i]  = datum -> z;
      }
      
      if (! dcycle (N, t(), r(), f(), period, area)) {
	cerr << "cycint: couldn't fit cycle " << ncycles << endl;
	return EXIT_FAILURE;
      }

      t1   = t[N - 2];  t2 = t[N - 1];
      r1   = r[N - 2];  r2 = r[N - 1];
      end  = t1 + (t1 - t2) * r1 / (r2 - r1);
      freq = 1.0 / period;

      favg = (favg * (ncycles - 1) + freq) / (double) ncycles;
      aavg = (aavg * (ncycles - 1) + area) / (double) ncycles;

      if (verbose) {
	sprintf (line, "%.4e %.4e %.4e %.4e %.4e",
		 end, freq, favg, area, aavg);
	cout << line << endl;
      }

      // -- Set up for next cycle.

      datum = new triplet (t1, r1, f1);
      data.push (datum);
      datum = new triplet (t2, r2, f2);
      data.push (datum);
      t1 = t2;
      r1 = r2;
      f1 = f2;
    }
  }  

  if (ncycles == 0) {
    cerr << "cycint: no second crossing" << endl;
    return EXIT_FAILURE;

  } else {
    sprintf (line, "averages:  frequency  %.4e,      area %.4e", favg, aavg);
    cout << line << endl;
  }

  return EXIT_SUCCESS;
}
  

static int dcycle (const int     np    ,
		   const double* x     ,
		   const double* y     ,
		   const double* z     ,
		   double&       period,
		   double&       area  )
// ---------------------------------------------------------------------------
// Given tables (length np) x (abscissa) and y (ordinate) of a function,
// estimate period and area of auxillary function z.
// Return 0 if no cycle can be found.
// ---------------------------------------------------------------------------
{
  double       xlo, xhi, xtp;
  const double EPS = 1.0e-6, HUGE = 1.0E99;

  // -- Initalize global variables.

  cc   = new double [np];
  xx   = (double*) x;
  yy   = (double*) y;
  num  = np;
  sign = 1.0;

  // -- Generate global natural spline coefficients cc.

  Recipes::spline (xx, yy, num, HUGE, HUGE, cc);

  // -- Find zero crossings.

  xlo = x[0];
  xhi = xlo + 0.05 * (x[np - 1] - x[0]);
  xlo = Recipes::rtsec (wave, xlo, xhi, EPS);

  xhi = x[np - 1];
  xtp = xhi - 0.05 * (xhi - xlo);
  xhi = Recipes::rtsec (wave, xtp, xhi, EPS);

  period = xhi - xlo;

  // -- Find area of auxillary function.

  yy = (double*) z;
  Recipes::spline (xx, yy, num, HUGE, HUGE, cc);

  area = dsplquad (x, z, cc, num, xlo, xhi);

  delete [] cc;

  return 1;
}


static double wave (const double& x)
// ---------------------------------------------------------------------------
// Return sign * interpolating spline function.
// ---------------------------------------------------------------------------
{
  double  y;

  Recipes::splint (xx, yy, cc, num, x, y);

  return sign * y;
}
