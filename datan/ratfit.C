///////////////////////////////////////////////////////////////////////////////
// ratfit.C: Chebyshev rational polynomial fit to data points, adaption
// of Numerical Recipes.
//
//               p0 + p1.x + p2.x^2 + ... pm.x^m
// f(x) \approx  -------------------------------
//               1  + q1.x + q2.x^2 + ... qk.x^k
//
// Usage: ratfit [options] [file]
// options:
//   -h      ... print this message
//   -n <mm> ... set order of numerator   to mm         [Default: 1]
//   -d <kk> ... set order of denominator to kk         [Default: 1]
//   -f <nn> ... print out nn evaluations from ll to hh [Default: 0]
//   -l <ll> ... lower limit for evaluations            [Default: 0]
//   -u <hh> ... upper limit for evaluations            [Default: 1]
//
// Input: x & y data points, at least kk + mm + 3 pairs.
// Output: coefficients of numerator and denominator polynomials.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cstdio>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

#include <Utility.h>
#include <Stack.h>
#include <Array.h>
#include "nr77.h"

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))

class doublet {
public:
  doublet (double Z, double W) : z (Z), w (W) { }
  double   z, w;
};

static char prog[]  = "ratfit";


void ratlsq (const double* x  ,
	     const double* y  ,
	     const int     npt,
	     double*       cof,
	     double&       dev,
	     const int     mm ,
	     const int     kk )
// ---------------------------------------------------------------------------
// Adaption of NR routine that accepts data as input.  Created matrices
// are transposed from NR format, in for use by FORTRAN subroutines.
// ---------------------------------------------------------------------------
{
  const int ncof = mm + kk + 1;
  if (npt < ncof) message ("ratlsq", "insufficient points", ERROR);

  char           buf[StrMax];
  const int      MAXIT = 5;
  register int   i, j, it;
  double         power, devmax, sum, e;  
  vector<double> work (3 * npt + 2 * ncof + (npt + ncof) * ncof);
  double*        bb   = work();
  double*        coff = bb   + npt;
  double*        ee   = coff + ncof;
  double*        u    = ee   + npt;
  double*        v    = u    + npt  * ncof;
  double*        w    = v    + ncof * ncof;
  double*        wt   = w    + ncof;

  for (i = 0; i < npt; i++) wt[i] = ee[i] = 1.0;
  dev = FLT_MAX;

  for (e = 0.0, it = 0; it < MAXIT; it++) {
    for (i = 0; i < npt; i++) {
      power = wt[i];
      bb[i] = power * (y[i] + SIGN (e, ee[i]));
      for (j = 0; j <= mm; j++) {
	u[j * npt + i] = power;
	power *= x[i];
      }
      power = -bb[i];
      for (j = mm + 1; j < ncof; j++) {
	power *= x[i];
	u[j * npt + i] = power;
      }
    }
    Recipes::svdcmp (u, npt, ncof, npt, ncof, w, v);
    Recipes::svbksb (u, w, v, npt, ncof, npt, ncof, bb, coff);
    devmax = sum = 0.0;
    for (j = 0; j < npt; j++) {
      ee[j] = Recipes::ratval (x[j], coff, mm, kk) - y[j];
      wt[j] = fabs (ee[j]);
      sum  += wt[j];
      if (wt[j] > devmax) devmax = wt[j];
    }
    e = sum / npt;
    if (devmax <= dev) {
      for (j = 0; j < ncof; j++) cof[j] = coff[j];
      dev = devmax;
    }
    
    sprintf (buf, ": iteration= %2d max error= %10.3e\n", it, devmax);
    message ("ratlsq", buf, REMARK);
  }
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  istream*        input;
  vector<double>  work;
  double          dev, xi, yi, *x, *y, *cof, *X, *Y;
  int             i, N, mm = 1, kk = 1, ncof;
  int             neval = 0;
  double          lo = 0.0, hi = 1.0;
  doublet*        datum;
  Stack<doublet*> data;

  while (--argc && **++argv == '-')
    switch (i = *++argv[0]) {
    case 'n':
      if (*++argv[0])
	mm = atoi (*argv);
      else {
	mm = atoi (*++argv);
	argc--;
      }
      break;
    case 'd':
      if (*++argv[0])
	kk = atoi (*argv);
      else {
	kk = atoi (*++argv);
	argc--;
      }
      break;
    case 'f':
      if (*++argv[0])
	neval = atoi (*argv);
      else {
	neval = atoi (*++argv);
	argc--;
      }
      break;
    case 'l':
      if (*++argv[0])
	lo = atof (*argv);
      else {
	lo = atof (*++argv);
	argc--;
      }
      break;
    case 'u':
      if (*++argv[0])
	hi = atof (*argv);
      else {
	hi = atof (*++argv);
	argc--;
      }
      break;
    case 'h': default:
      cout <<
	"Usage: ratfit [options] [file]\n"
	"options:\n"
	"  -h      ... print this message\n"
	"  -n <mm> ... set order of numerator   to mm         [Default: 1]\n"
	"  -d <kk> ... set order of denominator to kk         [Default: 1]\n"
	"  -f <nn> ... print out nn evaluations from ll to hh [Default: 0]\n"
	"  -l <ll> ... lower limit for evaluations            [Default: 0]\n"
	"  -u <hh> ... upper limit for evaluations            [Default: 1]\n";
      return EXIT_SUCCESS;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) message (prog, "unable to open input file", ERROR);
  } else input = &cin;

  while (*input >> xi >> yi) {
    datum = new doublet (xi, yi);
    data.push (datum);
  }

  work.setSize (2 * (N = data.depth()) + (ncof = mm + kk + 1) + 2 * neval);
  x   = work();
  y   = x   + N;
  cof = y   + N;
  X   = cof + ncof;
  Y   = X   + neval;

  for (i = 0; i < N; i++) {
    datum = data.pop();
    x[N-i-1] = datum -> z;
    y[N-i-1] = datum -> w;
  }

  ratlsq (x, y, N, cof, dev, mm, kk);

  cout << "Numerator:" << endl;
  for (i = 0; i <= mm; i++) cout << cof[i] << endl;

  cout << "Denominator:" << endl;

  cout << "1.0" << endl;

  for (i = mm + 1; i < ncof; i++) cout << cof[i] << endl;

  if (neval)
    for (i = 0; i < neval; i++) {
      X[i] = lo + (hi - lo) * i/(neval - 1);
      Y[i] = Recipes::ratval (X[i], cof, mm, kk);
      cout << X[i] << "\t\t" << Y[i] << endl;
    }

  return EXIT_SUCCESS;
}
