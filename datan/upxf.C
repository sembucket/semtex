///////////////////////////////////////////////////////////////////////////////
// upxf.C: compute zero-upcrossing frequency, and peak values of oscillatory
// waveform.
//
// Usage:
// ------
// upxf [-h] [-v] [file]
//
// Synopsis:
// ---------
// Expected input is 2 columns of ASCII, the first is time, the second is
// some process.  Compute the average upcrossing frequency, print up.
// Find zeros by iteration, maxima/minima by optimization techniques, both
// based on cubic spline fits to data cycles.
//
// Options:
// --------
// -h (help):    print usage string.
// -v (verbose): peaks, local and running average frequencies.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <strstream.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>

#include <Utility.h>
#include <Stack.h>
#include <Array.h>
#include <nr77.h>

class doublet {
public:
  doublet (double X, double Y) : x (X), y (Y) { }
  double    x, y;
};

static int     num;
static double  sign;
static double* xx;
static double* yy;
static double* cc;
static double  wave   (const double&);
static int     dcycle (const int, const double*, const double*, 
		       double&, double&, double&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  ifstream        file;
  doublet*        datum;
  Stack<doublet*> data;
  vector<double>  t, f;
  int             i, N;

  double t1, t2, f1, f2, period, end;
  double freq, favg = 0.0;
  double fmax, fmav = 0.0;
  double fmin, fmiv = 0.0;

  char   line[StrMax], c;
  int    crossed = 0, verbose = 0, ncycles = 0;
  char   usage[] = "upxf [-h] [-r] [-v] [file]";

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

  if   (argc == 1) file.open   (*argv, ios::in);
  else             file.attach (0);

  if (!file) {
    cerr << "unable to open input file" << endl;
    return EXIT_FAILURE;
  }

  // -- Find first crossing.

  file >> t1 >> f1;
  while (!crossed && file >> t2 >> f2) {
    if (crossed = ((f1 < 0.0) && (f2 >= 0.0))) {
      datum = new doublet (t1, f1);
      data.push (datum);
      if (verbose) {
	end = t1 + (t1 - t2) * f1 / (f2 - f1);
	sprintf (line, "%.4e", end);
	cout << line << endl;
      }
    } 
    t1 = t2;
    f1 = f2;
  }

  if (!crossed) {
    cerr << "upxf: no crossing found" << endl;
    return EXIT_FAILURE;

  } else {
    datum = new doublet (t2, f2);
    data.push (datum);
    crossed = 0;
  }

  // -- Do the rest of input, cycle by cycle.

  while (file) {

    while (!crossed && file >> t2 >> f2) {
      crossed = ((f1 < 0.0) && (f2 >= 0.0));
      datum = new doublet (t2, f2);
      data.push (datum);
      t1 = t2;
      f1 = f2;
    }

    if (crossed) {		// -- End of current cycle.
      i = N = data.depth();
      t.setSize (N);
      f.setSize (N);

      ncycles += 1;
      crossed  = 0;

      while (i--) {
	datum = data.pop();
	t[i]  = datum -> x;
	f[i]  = datum -> y;
      }
      
      if (! dcycle (N, t(), f(), period, fmin, fmax)) {
	cerr << "upxf: couldn't fit cycle " << ncycles << endl;
	return EXIT_FAILURE;
      }

      t1   = t[N - 2];  t2 = t[N - 1];
      f1   = f[N - 2];  f2 = f[N - 1];
      end  = t1 + (t1 - t2) * f1 / (f2 - f1);
      freq = 1.0 / period;

      favg = (favg * (ncycles - 1) + freq) / (double) ncycles;
      fmav = (fmav * (ncycles - 1) + fmax) / (double) ncycles;
      fmiv = (fmiv * (ncycles - 1) + fmin) / (double) ncycles;

      if (verbose) {
	sprintf (line, "%.4e %.4e %.4e %.4e %.4e %.4e %.4e",
		 end, freq, favg, fmin, fmiv, fmax, fmav);
	cout << line << endl;
      }

      // -- Set up for next cycle.

      datum = new doublet (t1, f1);
      data.push (datum);
      datum = new doublet (t2, f2);
      data.push (datum);
      t1 = t2;
      f1 = f2;
    }
  }  

  if (ncycles == 0) {
    cerr << "upxf: no second crossing" << endl;
    return EXIT_FAILURE;

  } else {
    sprintf (line, 
	     "averages:  frequency  %.4e,        min %.4e,       max %.4e",
	     favg, fmiv, fmav);
    cout << line << endl;
  }

  return EXIT_SUCCESS;
}
  

static int dcycle (const int     np    ,
		   const double* x     ,
		   const double* y     ,
		   double&       period,
		   double&       min   ,
		   double&       max   )
// ---------------------------------------------------------------------------
// Given tables (length np) x (abscissa) and y (ordinate) of a function,
// estimate extrema and period.  Return 0 if no cycle can be found.
// ---------------------------------------------------------------------------
{
  double       xlo, xhi, xav, xtp;
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

  // -- Find turning points.

  xav = 0.5 * (xlo + xhi);
  xtp = 0.5 * (xlo + xav);
  
  if (wave (xtp) > 0.0) sign = -1.0;

  min = sign * Recipes::golden (xlo, xtp, xav, wave, EPS, xtp);

  xtp  = 0.5 * (xav + xhi);
  sign = -sign;
  max = sign * Recipes::golden (xav, xtp, xhi, wave, EPS, xtp);

  if (sign > 0.0) {
    double tmp = min;
    min = max;
    max = tmp;
  }

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
