///////////////////////////////////////////////////////////////////////////////
// momenta.C: from a file of data (optionally stdin) that lists
// displacement and weights, compute the first four moments of the
// weights.
// 
// Displacements and weights are taken to be drawn from
// spatially-periodic data.
//
// Usage: momenta [file]
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <strstream.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include <Utility.h>
#include <Stack.h>
#include <Array.h>

class doublet {
public:
  doublet (double Z, double W) : z (Z), w (W) { }
  double   z, w;
};

static char prog[]  = "momenta";


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  double          x, y, zpk, zav, wt, ww, dz, lz, shift;
  double          sum, nfac, mean, sdev, var, skew, flat;
  double          wmax = -FLT_MAX, zmin = FLT_MAX, zmax = -FLT_MAX;
  doublet*        datum;
  Stack<doublet*> data;
  vector<double>  z, w;
  int             i, N;

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h': default:
      cout << prog << " [file]" << endl;
      return EXIT_SUCCESS;
    }
    
  if (argc == 1) {
    ifstream* inputfile = new ifstream (*argv);
    if (inputfile -> good()) {
      cin = *inputfile;
      } else {
	cerr << prog << ": unable to open input file" << endl;
	exit (EXIT_FAILURE);
    }
  }

  while (cin >> x >> y) {
    zmin  = min (x, zmin);
    zmax  = max (x, zmax);
    datum = new doublet (x, y);
    data.push (datum);
  }

  z.setSize (N = data.depth());
  w.setSize (N);

  lz  = zmax - zmin;
  dz  = lz / (N - 1.0);
  lz += dz;

  for (sum = 0.0, i = 0; i < N; i++) {
    datum = data.pop();
    x = datum -> z;
    y = datum -> w;
    if (fabs (y) > wmax) { wmax = fabs (y); zpk  = x; }
    sum     += y;
    z[N-i-1] = x;
    w[N-i-1] = y;
  }

  nfac = 1.0 / sum;

  // -- Shift z locations so that peak is roughly centered in domain.

  shift = zpk - 0.5*lz;
  for (i = 0; i < N; i++) z[i] = fmod (z[i] - shift + lz, lz);

  // -- Find z location of mean.

  for (zav = 0.0, i = 0; i < N; i++) zav += z[i] * w[i];
  zav *= nfac;

  // -- Compute higher moments.

  for (var = 0.0, skew = 0.0, flat = 0.0, i = 0; i < N; i++) {
    wt    = z[i] - zav;
    ww    = wt * wt * w[i];
    var  += ww;
    ww   *= wt;
    skew += ww;
    ww   *= wt;
    flat += ww;
  }
  
  mean  = zav + shift;
  var  *= nfac;
  skew *= nfac;
  flat *= nfac;
  
  sdev  = sqrt (var);
  skew  = skew / (var * sdev);
  flat  = flat / (var * var ) - 3.0;

  cout.precision (8);
  cout
    << mean << "  "
    << sdev << "  "
    << var  << "  "
    << skew << "  "
    << flat << endl;

  return EXIT_SUCCESS;
}
