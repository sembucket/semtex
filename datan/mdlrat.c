/*****************************************************************************
 * mdlrat.c: produce ensemble-averaged value of the ratio of sums
 * of modal energies normalized by zero mode's energy.  Input is a 
 * single column of modal energies, the block size is nominated
 * on the command line.
 *
 * Usage:
 * ------
 * mdlrat [-h] [-b <#block>] [-v] [input]
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char prog[] = "mdlrat";

static void getargs    (int, char**, int*, int*, FILE**);
static int  readdata   (FILE*, const int, double*);
static void energyrat  (const int, const double*, double*, double*);

#define DEFAULT 72


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  FILE   *fp_in  = stdin,
         *fp_out = stdout;
  int    verb = 0, nblk = DEFAULT, navg = 0;
  double zero = 0.0, sum = 0.0, ratio = 0.0; 
  double *data;


  /* -- Parse command line. */

  getargs (argc, argv, &nblk, &verb, &fp_in);

  /* -- Set up storage. */

  data = (double*) calloc (nblk, sizeof (double));
  
  /* -- Accumulate autocorrelations. */

  while (readdata (fp_in, nblk, data) == nblk) {
    energyrat   (nblk, data, &zero, &sum);
    ratio += sum / zero;
    navg++;
  }
     
  /* -- Scale, print, exit. */

  printf ("%g\n", ratio / navg);

  return EXIT_SUCCESS;
}


static void getargs (int    argc  ,
		     char** argv  ,
		     int*   nblk  ,
		     int*   verb  ,
		     FILE** fp_in )
/* ------------------------------------------------------------------------- *
 * Parse command line.
 * ------------------------------------------------------------------------- */
{
  char   c;
  char   line[FILENAME_MAX];
  static char usage[] = 
    "usage: %s [options] [input]\n"
    "[options]:\n"
    "-h          ... display this message\n"
    "-b <#block> ... size of data block         [Default: %1d]\n"
    "-v          ... toggle verbose mode\n";

  while (--argc && (*++argv)[0] == '-') {
    switch (c = *++argv[0]) {
    case 'h':
      fprintf (stderr, usage, prog, DEFAULT);
      exit    (EXIT_SUCCESS);
      break;
    case 'b':
      if (*++argv[0])
	*nblk = atoi (*argv);
      else {
        *nblk = atoi (*++argv);
	argc--;
      }
      break;
    case 'v':
      *verb = 1;
      break;
    default:
      fprintf (stderr, "%s: '%c': unknown option\n", prog, c);
      fprintf (stderr, usage, prog, DEFAULT);
      exit    (EXIT_FAILURE);
      break;
    }
  }

  if (argc == 1 && !(*fp_in = fopen (*argv, "r"))) {
    fprintf (stderr, "%s: couldn't open file %s\n", prog, *argv);
    exit    (EXIT_FAILURE);
  }
}


static int readdata (FILE*     fp  ,
		     const int nblk,
		     double*   data)
/* ------------------------------------------------------------------------- *
 * Try to read block of data length nblk into data, return number read.
 * ------------------------------------------------------------------------- */
{
  register int i;

  for (i = 0; i < nblk; i++) if (fscanf (fp, "%lf", data + i) != 1) break;

  return i;
}


static void energyrat (const int     np,
		       const double* x ,
		       double*       mz,
		       double*       ms)
/* ------------------------------------------------------------------------- *
 * Compute in c the circular correlation of x & y.  X & y are np
 * (number of points) long, while c is at least nl (number of lags)
 * long.  Nl can be greater than np.  Mean values of x & y are
 * returned in *mx & *my.
 * ------------------------------------------------------------------------- */
{
  register int    i;
  register double s = 0.0;

  for (i = 1; i < np; i++) s += x[i];
  
  *mz = x[0];
  *ms = s;
}
