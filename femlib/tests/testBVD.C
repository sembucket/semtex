//////////////////////////////////////////////////////////////////////////////
// filter.C: Generate/test filtering functions.
//
// Usage: filter -a <num> -N <num> -p <num> -s <num>
// where
// -N <num> ... supplies the number of points in the filter [0, N]
// -s <num> ... supplies the filter lag
// -p <num> ... supplies the filter order
// -a <num> ... supplies attenuation factor at high frequencies
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


static void bvdFilter (const int N,
		       const int p,
		       const int s,
		       real      a,
		       real*     F)
// ---------------------------------------------------------------------------
// Load F with the Boyd--Vandeven filter [0, N] of order p, lag s.
// NB: N should be one less than the number of coefficients to
// which the filter will be applied.
// 
// Input parameter t gives the attenuation at high wavenumbers. 0<=t<=1,
// with t = 1 giving complete attenuation.
// ---------------------------------------------------------------------------
{
  int        i;
  real       arg, theta, chi, omega;
  const real EPS = EPSSP;

  for (i = 0; i < s; i++)
    F[i] = 1.0;
  
  for (i = s; i <= N; i++) {
    theta = (i - s) / (real) (N - s);
    omega = fabs(theta) - 0.5;
    if ((fabs (theta - 0.5)) < EPS) 
      chi = 1.0;
    else {
      arg = 1.0 - 4.0 * sqr (omega);
      chi = sqrt (-log (arg) / (4.0 * sqr (omega)));
    }
    F[i] = (1.0 - a) + a * 0.5 * erfc (2.0*sqrt(p)*chi*omega);
  }
}

static void getargs (int     argc,
		     char**  argv,
		     real&   a   ,
		     int&    s   ,
		     int&    p   ,
		     int&    N   )
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  const char *usage = 
    "Usage: filter -a <num> -N <num> -p <num> -s <num>\n"
    "where\n"
    "-a <num> ... supplies attenuation factor at high frequencies [0, 1]\n"
    "-N <num> ... supplies the number of points in the filter [0, N]\n"
    "-s <num> ... supplies the filter lag\n"
    "-p <num> ... supplies the filter order\n";

  if (argc != 9) {
    cerr << usage;
    exit (EXIT_FAILURE);
  }
  
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'a':
      --argc;
      a = atof (*++argv);
      break;
    case 'N':
      --argc;
      N = atoi (*++argv);
      break;
    case 'p':
      --argc;
      p = atoi (*++argv);
      break;
    case 's':
      --argc;
      s = atoi (*++argv);
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
  int  i, s, p, N;
  real a;

  getargs (argc, argv, a, s, p, N);

  vector<real> work (N + 1);
  real         *filter = work();

  bvdFilter (N, p, s, a, filter);

  for (i = 0; i <= N; i++) 
    cout << i << '\t' << filter[i] << endl;
  
  return EXIT_SUCCESS;
}
