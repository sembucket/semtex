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
  ifstream       file;
  double         x, y, zmax, wmax;
  doublet        datum;
  Stack<doublet> data;
  vector<double> ztmp, z, wtmp, w;
  int            i, ihalf, imax, N, shift;

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h': default:
      cout << prog << " [file]" << endl;
      return EXIT_SUCCESS;
      break;
    }

  if   (argc == 1) file.open   (*argv, ios::in);
  else             file.attach (0);

  if (!file) {
    cerr << prog << ": unable to open input file" << endl;
    return EXIT_FAILURE;
  }

  while (file >> x >> y) {
    datum = new doublet (x, y);
    data.push (datum);
  }

  z   .setSize (N = data.depth());
  ztmp.setSize (N);
  w   .setSize (N);
  wtmp.setSize (N);

  for (i = 0; i < N; i++) {
    datum = data.pop();
    ztmp[N-i-1] = datum -> z;
    wtmp[N-i-1] = datum -> w;
  }

  // -- Rearrange data to place maximum at central location in storage.

  for (wmax = -FLT_MAX, i = 0; i < N; i++)
    if (fabs (wtmp[i]) > wmax) {
      wmax = fabs (wtmp[i]);
      zmax = ztmp[i];
      imax = i;
    }

  shift = imax - N / 2;

  for (i = 0; i < N; i++) {
    z [(i + N - shift) % N] = ztmp[i];
    w [(i + N - shift) % N] = wtmp[i];
  }

  
  return EXIT_SUCCESS;
}
