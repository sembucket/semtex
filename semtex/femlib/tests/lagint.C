//////////////////////////////////////////////////////////////////////////////
// lagint.C: given a table of N values, assumed located at the GLL
// points in [-1, 1], return I values of the N-1 polynomial
// interpolant for locations uniformly distributed on the interval,
// computed using the Lagrange polynomials.
//
// Usage: lagint -i <num> [file]
// where
// -i <num> ... number of interpolant points on [-1, 1]
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>

#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>

#include <femdef.h>
#include <Array.h>
#include <Veclib.h>
#include <Femlib.h>
#include <Blas.h>
#include <Utility.h>
#include <Stack.h>


void getargs (int    argc,
	      char** argv,
	      int&   nint)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: lagint -i <num> [file]\n"
    "where\n"
    "-i <num> ... number of interpolant points on [-1, 1]\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'i':
      if (*++argv[0]) nint = atoi (*argv);
      else { --argc;  nint = atoi (*++argv); }
      break;
    default:
      cerr << usage; exit (EXIT_FAILURE);
      break;
    }

  if (!nint) { cerr << usage; exit (EXIT_FAILURE); }

  if (argc == 1) {
    ifstream* inputfile = new ifstream (*argv);
    if (inputfile -> good()) {
      cin = *inputfile;
      } else {
	cerr <<  "lagint: unable to open input file" << endl;
	exit (EXIT_FAILURE);
    }
  }
}


int loadVals (istream&        file,
	      vector<double>& val )
// ---------------------------------------------------------------------------
// Get nodal values from file, return number of values.
// ---------------------------------------------------------------------------
{
  int           ntot, num = 0;
  double        datum;
  Stack<double> data;

  while (file >> datum) { data.push (datum); num++; }

  val.setSize (ntot = num);

  while (num--) val[num] = data.pop();

  return ntot;
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int            i, j, ngll, nint = 0;
  vector<double> u, v, w, x, z;
  double         **II, **IT;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, nint);

  ngll = loadVals (cin, v);
  z.setSize (ngll);
  w.setSize (ngll);

  x.setSize (nint);
  u.setSize (nint);

  for (i = 0; i < nint; i++) x[i] = -1.0 + i * 2.0/(nint - 1);

  II = dmatrix (0, nint-1, 0, ngll-1);
  IT = dmatrix (0, ngll-1, 0, nint-1);

  Femlib::GLLzw       (ngll, z(), w());
  Femlib::LagrangeInt (ngll, z(), nint, x(), II, IT);

  Blas::mxv (*II, nint, v(), ngll, u());

  for (i = 0; i < nint; i++) cout << x[i] << '\t' << u[i] << endl;

  return EXIT_SUCCESS;
}
