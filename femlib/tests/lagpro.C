//////////////////////////////////////////////////////////////////////////////
// lagpro.C: given a table of N values, assumed located at the GLL
// points in [-1, 1], return values of the original interpolating
// polynomial, order N-1, at I GLL points; also return the values of
// the order I-1 interpolating polynomial for the new I GLL points
// evaluated at the original N GLL points.
//
// Usage: lagpro -i <num> [file]
// where
// -i <num> ... number of new GLL points on [-1, 1]
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


void getargs (int       argc,
	      char**    argv,
	      int&      nint,
	      ifstream& file)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: lagpro -i <num> [file]\n"
    "where\n"
    "-i <num> ... number of new GLL points on [-1, 1]\n";

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
    
  if   (argc == 1) file.open   (*argv, ios::in);
  else             file.attach (0);

  if (!file) {
    cerr << "lagint: unable to open input file" << endl; exit (EXIT_FAILURE);
  }
}


int loadVals (ifstream&       file,
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
  ifstream       file;
  int            i, j, nold, nnew = 0;
  vector<double> u, v;
  const double   **IF, **IB;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, nnew, file);

  nold = loadVals (file, v);

  u.setSize (nnew);

  if (nold == nnew)
    Veclib::copy (nold, v(), 1, u(), 1);
  else {
    Femlib::mesh (GLL, GLL, nold, nnew, 0, &IF, 0, 0, 0);
    Femlib::mesh (GLL, GLL, nnew, nold, 0, &IB, 0, 0, 0);
    Blas::mxv    (*IF, nnew, v(), nold, u());
    Blas::mxv    (*IB, nold, u(), nnew, v());
  }

  for (i = 0; i < nnew; i++) cout << i << '\t' << u[i] << endl;
  cout << endl;
  for (i = 0; i < nold; i++) cout << i << '\t' << v[i] << endl;

  return EXIT_SUCCESS;
}
