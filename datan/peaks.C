///////////////////////////////////////////////////////////////////////////////
// peaks.C: find and print minima or maxima of a function, input as a
// 2-column file: first column is a time/position value, second is
// function value.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdlib.h>


int main (int    argc,
	  char** argv)
{
  ifstream file;
  double   func[3], time[3];
  int      min = 0, verbose = 0;
  char usage[] = "peaks [-h] [-min] [-v] [file]", c;

  // -- Process command line.

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      cout << usage << endl;
      return EXIT_SUCCESS;
      break;
    case 'm':
      if (argv[0][1] == 'i') min = 1;
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

  // -- Initialise data windows.

  file >> time[2] >> func[2];
  file >> time[1] >> func[1];
  if (!file)  {
    cerr << "unable to initialise from input file" << endl;
    return EXIT_FAILURE;
  }
  
  cout << setprecision(8);

  // -- Main loop.

  while (file >> time[0] >> func[0]) {
    if (min) {
      if (func[1] <  func[0] && func[1] <  func[2] ||
	  func[1] <  func[0] && func[1] == func[2])
      cout << time[1] << "  " << func[1] << endl;
    } else {
      if (func[1] >  func[0] && func[1] >  func[2] ||
	  func[1] >  func[0] && func[1] == func[2])
      cout << time[1] << '\t' << func[1] << endl;
    }
    time[2] = time[1]; func[2] = func[1];
    time[1] = time[0]; func[1] = func[0];
  }
  
  return EXIT_SUCCESS;
}

