//////////////////////////////////////////////////////////////////////////////
// glzw.C: print out Gauss--Lobatto nodes and weights on [-1, 1].
//
// For now, just for the Legendre basis.
//
// Usage: glzw -N <num>
// where
// -N <num> ... supplies the number of GL nodes.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>

#include <iostream.h>
#include <iomanip.h>

#include <femdef.h>
#include <Array.h>
#include <Veclib.h>
#include <Femlib.h>
#include <Blas.h>
#include <Utility.h>

static void getargs (int     argc,
		     char**  argv,
		     int&    N   )
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: glzw -N <num>\n"
    "where\n"
    "-N <num> ... supplies the number of GL nodes\n";

  if (argc != 3) {
    cerr << usage;
    exit (EXIT_FAILURE);
  }
  
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'N':
      if (*++argv[0]) N = atoi (*argv);
      else { --argc;  N = atoi (*++argv); }
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int i, N;

  getargs (argc, argv, N);

  vector<real> work (2*N);
  real         *z = work();
  real         *w = z + N;

  Femlib::GLLzw (N, z, w);

  cout.precision (8);
  cout.setf (ios::fixed, ios::floatfield);
  for (i = 0; i < N; i++) 
    cout << i << '\t' << z[i] << '\t' << w[i] << endl;
  
  return EXIT_SUCCESS;
}
