///////////////////////////////////////////////////////////////////////////////
// cycle.C: print out evenly spaced samples of a single period of a waveform.
//
// Usage:
// ------
// cycle [-h] [-n num] [file]
//
// Synopsis:
// ---------
// Expected input is 3 columns of ASCII, the first is time, the second is
// a process whose zero crossings are used to define cycle times of the 
// process in the third column.  Find the first cycle of the process in the
// second column, interpolate equally-spaced values of all three columns
// using cubic splines, and print up the results.
//
// Options:
// --------
// -h    : Print usage string.
// -n num: Interpolate num equally-spaced values (Default: 512).
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <strstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cfloat>

using namespace std;

#include <Utility.h>
#include <Veclib.h>
#include <Stack.h>
#include <Array.h>
#include "nr77.h"

class triplet {
public:
  triplet (double X, double Y, double Z) : x (X), y (Y), z (Z) { }
  double    x, y, z;
};

static const char prog[]  = "cycle";
static const int  DEFAULT = 512;
static int     num;
static double* xx;
static double* yy;
static double* cc;
static double  wave   (const double&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  istream*        input;
  triplet*        datum;
  Stack<triplet*> data;
  vector<double>  t, r, f;	// -- Time, reference, function.
  vector<double>  rcof, fcof;	// -- Spline coefficients for r & f.
  int             i, N, nsamp = DEFAULT;
  double          t1, t2, r1, r2, f1, f2, tlo, thi, period;
  int             crossed = 0, ncycles = 0;
  char            usage[] = ": [-h] [-n num] [file]";

  // -- Process command line.

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << prog << usage << endl;
      return EXIT_SUCCESS;
      break;
    case 'n':
      if (*++argv[0])
	nsamp = atoi (*argv);
	else {
	  nsamp = atoi (*++argv);
	  argc--;
	}
      break;
    default:
      cerr << usage << endl;
      return EXIT_FAILURE;
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) message (prog, "unable to open input file", ERROR);
  } else input = &cin;

  // -- Find first crossing.

  *input >> t1 >> r1 >> f1;
  while (!crossed && *input >> t2 >> r2 >> f2) {
    if (crossed = ((r1 < 0.0) && (r2 >= 0.0))) {
      datum = new triplet (t1, r1, f1);
      data.push (datum);
    } 
    t1 = t2;
    r1 = r2;
    f1 = f2;
  }

  if (!crossed) {
    cerr << prog << ": no crossing found" << endl;
    return EXIT_FAILURE;

  } else {
    datum = new triplet (t2, r2, f2);
    data.push (datum);
    crossed = 0;
  }

  // -- Do first cycle of input found.

  while (!crossed && *input >> t2 >> r2 >> f2) {
    crossed = ((r1 < 0.0) && (r2 >= 0.0));
    datum = new triplet (t2, r2, f2);
    data.push (datum);
    t1 = t2;
    r1 = r2;
    f1 = f2;
  }

  // -- Should have first cycle.

  if (crossed) {
    i = N = data.depth();
    t.setSize    (N);
    r.setSize    (N);
    f.setSize    (N);
    rcof.setSize (N);
    fcof.setSize (N);

    ncycles += 1;
    crossed  = 0;

    while (i--) {
      datum = data.pop();
      t[i]  = datum -> x;
      r[i]  = datum -> y;
      f[i]  = datum -> z;
    }

    // -- Generate natural spline coefficients for r & f.

    Veclib::spline (N, FLT_MAX, FLT_MAX, t(), r(), rcof());
    Veclib::spline (N, FLT_MAX, FLT_MAX, t(), f(), fcof());

    // -- Set up for root-finding.

    num = N;
    xx  = (double*) t();
    yy  = (double*) r();
    cc  = (double*) rcof();

    // -- Find zero crossings of r.

    tlo = t[0];
    thi = tlo + 0.05 * (t[N - 1] - t[0]);
    tlo = Recipes::rtsec (wave, tlo, thi, EPSm6);

    thi = t[N - 1];
    t2  = thi - 0.05 * (thi - tlo);
    thi = Recipes::rtsec (wave, t2,  thi, EPSm6);

    period = thi - tlo;
    t2     = period / nsamp;
    
    for (i = 0; i < nsamp; i++) {
      t1 = tlo + i * t2;
      r1 = Veclib::splint (N, t1, t(), r(), rcof());
      f1 = Veclib::splint (N, t1, t(), f(), fcof());
      cout << setw(16) << t1 << setw(16) << r1 << setw(16) << f1 << endl;
    }
  }  

  if (ncycles == 0) {
    cerr << prog << ": no second crossing" << endl;
    return EXIT_FAILURE;

  }

  return EXIT_SUCCESS;
}


static double wave (const double& x)
// ---------------------------------------------------------------------------
// Return interpolating spline function evaluated at x.
// ---------------------------------------------------------------------------
{
  double y;

  Recipes::splint (xx, yy, cc, num, x, y);
  return y;
}
