//////////////////////////////////////////////////////////////////////////////
// lagder.C: given a table of N values, assumed located at the GLL
// points in [-1, 1], return values of the derivative at the same points.
//
// Usage: lagder [file]
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
	      ifstream& file)
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: lagder [file]\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    default:
      cerr << usage; exit (EXIT_FAILURE);
      break;
    }

  if   (argc == 1) file.open   (*argv, ios::in);
  else             file.attach (0);

  if (!file) {
    cerr << "lagder: unable to open input file" << endl; exit (EXIT_FAILURE);
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
  int            i, np;
  vector<double> u, v;
  const double   **DV;

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);

  getargs (argc, argv, file);

  np = loadVals (file, v);

  u.setSize (np);

  Femlib::mesh (GLL, GLL, np, np, 0, 0, 0, &DV, 0);
  Blas::mxv    (*DV, np, v(), np, u());

  for (i = 0; i < np; i++) cout << i << '\t' << u[i] << endl;

  return EXIT_SUCCESS;
}
