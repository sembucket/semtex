///////////////////////////////////////////////////////////////////////////////
// zeros.C: find and print zero crossings of a function, input as a
// 2-column file: first column is a time/position value, second is
// function value.  Linear interpolation.
//
// Optionally check for crossings at level nominated on command line.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

static char prog[] = "zeros";


int main (int    argc,
	  char** argv)
{
  istream* input;
  double   x[2], y[2], zero, val = 0.0;
  int      verbose = 0;
  char     usage[] = "zeros [-h] [-v] [-z <num>] [file]";

  // -- Process command line.

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage << endl;
      return EXIT_SUCCESS;
      break;
    case 'v':
      verbose = 1;
      break;
    case 'z':
      if   (*++argv[0])           val = atof (*argv);
      else              { --argc; val = atof (*++argv); }
      break;
    default:
      cerr << usage << endl;
      return EXIT_FAILURE;
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) {
      cerr << "peaks: unable to open file" << endl; return EXIT_FAILURE;
    }
  } else input = &cin;

  // -- Initialise data windows.

  *input >> x[0] >> y[0];
  if (!*input) {
    cerr << "unable to initialise from input file" << endl;
    return EXIT_FAILURE;
  }
  y[1] -= val;
  
  cout << setprecision(8);

  // -- Main loop.

  while (*input >> x[1] >> y[1]) {
    y[1] -= val;
    zero  = x[0] + (x[1]-x[0])*y[0]/(y[0]-y[1]);
    if        (y[0] < 0.0 && y[1] >= 0.0) { // -- +ve xing.
      if (verbose) cout << "[+] ";
      cout << zero << endl;
    } else if (y[0] > 0.0 && y[1] <= 0.0) { // -- -ve xing.
      if (verbose) cout << "[-] ";
      cout << zero << endl;
    }
    x[0] = x[1]; y[0] = y[1];
  }
  
  return EXIT_SUCCESS;
}

