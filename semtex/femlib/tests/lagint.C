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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

#include <cfemdef>
#include <Array.h>
#include <veclib_h>
#include <femlib_h>
#include <blas_h>
#include <utility_h>
#include <Stack.h>

static char prog[] = "lagint";

void getargs (int       argc,
	      char**    argv,
	      int&      nint,
	      istream*& file)
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
    file = new ifstream (*argv);
    if (file -> bad()) message (prog, "unable to open input file", ERROR);
  } else file = &cin;
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
  istream*       input;
  int            i, j, ngll, nint = 0;
  vector<double> u, v, w, x, z;
  double         **II, **IT;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, nint, input);

  ngll = loadVals (*input, v);
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
