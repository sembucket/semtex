/*****************************************************************************
 * pcor.c: produce ensemble-averaged circular autocorrelation
 * function from columnar format input.
 *
 * Synopsis:
 * ---------
 * pcor reads ASCII data from stdin or named file. The data is read
 * in blocks of given size nominated on command line.  For each of
 * these blocks, the circular convolution is computed over a nominated
 * number of lags.  The circular correlations and the block means are
 * accumulated.  After all data are read, the product of the accumulated 
 * block means is subtracted from the accumulated circular correlations
 * and the outcome is normalised so the first value is unity.
 *
 * Usage:
 * ------
 * pcor [-h] [-b <#block>] [-n <#lags>] [-o outfile] [-v] [input]
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char prog[] = "pcor";

static void getargs    (int, char**, int*, int*, int*, FILE**, FILE**);
static int  readdata   (FILE*, const int, double*);
static void pcor       (const int, const int, const double*, const double*,
			double*, double*, double*);
static void accumulate (const int, const double*, double*,
			const double, double*);
static void normalize  (const int, const int, double*, const double);
static void printup    (FILE*, const int, const double*);

#define DEFAULT 32


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  FILE   *fp_in  = stdin,
         *fp_out = stdout;
  int    verb = 0, nblk = 2 * DEFAULT, nlag = DEFAULT, navg = 0;
  double mean, msum = 0.0;
  double *data, *csum, *acor;


  /* -- Parse command line. */

  getargs (argc, argv, &nlag, &nblk, &verb, &fp_in, &fp_out);

  /* -- Set up storage. */

  data = (double*) calloc (nblk + nlag + nlag, sizeof (double));
  acor = data + nblk;
  csum = acor + nlag;
  
  /* -- Accumulate autocorrelations. */

  while (readdata (fp_in, nblk, data) == nblk) {
    pcor       (nblk, nlag, data, data, acor, &mean, &mean);
    accumulate (nlag, acor, csum, mean, &msum);
    navg++;
  }
     
  /* -- Scale, print, exit. */

  normalize (nlag, navg, csum, msum);
    
  printup (fp_out, nlag, csum);

  return EXIT_SUCCESS;
}


static void getargs (int    argc  ,
		     char** argv  ,
		     int*   nlag  ,
		     int*   nblk  ,
		     int*   verb  ,
		     FILE** fp_in ,
		     FILE** fp_out)
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
    "-n <#lags>  ... compute this many lags     [Default: %1d]\n"
    "-o <output> ... write output to named file [Default: stdout]\n"
    "-v          ... toggle verbose mode\n";

  while (--argc && (*++argv)[0] == '-') {
    switch (c = *++argv[0]) {
    case 'h':
      fprintf (stderr, usage, prog, 2 * DEFAULT, DEFAULT);
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
    case 'n':
      if (*++argv[0])
	*nlag = atoi (*argv);
      else {
        *nlag = atoi (*++argv);
	argc--;
      }
      break;
    case 'o':
      if (*++argv[0])
        *fp_out = fopen (*argv, "w");
      else {
        *fp_out = fopen (*++argv, "w");
        argc--;
      }
      break;
    case 'v':
      *verb = 1;
      break;
    default:
      fprintf (stderr, "%s: '%c': unknown option\n", prog, c);
      fprintf (stderr, usage, prog, 2 * DEFAULT, DEFAULT);
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


static void pcor (const int     np,
		  const int     nl,
		  const double* x ,
		  const double* y ,
		  double*       c ,
		  double*       mx,
		  double*       my)
/* ------------------------------------------------------------------------- *
 * Compute in c the circular correlation of x & y.  X & y are np
 * (number of points) long, while c is at least nl (number of lags)
 * long.  Nl can be greater than np.  Mean values of x & y are
 * returned in *mx & *my.
 * ------------------------------------------------------------------------- */
{
  register int    i, j, l;
  register double ax = 0.0, ay = 0.0, corr;

  for (i = 0; i < np; i++) {
    ax += x[i];
    ay += y[i];
  }
  *mx = ax / np;
  *my = ay / np;

  for (l = 0; l < nl; l++) {
    corr = 0.0;
    for (i = 0; i < np; i++) {
      j = (i + l) % np;
      corr += x[i] * y[j];
    }
    c[l] = corr / np;
  }
}


static void accumulate (const int     nlag,
			const double* acor,
			double*       csum,
			const double  mean,
			double*       msum)
/* ------------------------------------------------------------------------- *
 * Accumulate buffer.
 * ------------------------------------------------------------------------- */
{
  register int i;

  for (i = 0; i < nlag; i++) csum[i] += acor[i];

  *msum += mean;
}


static void normalize  (const int    nlag,
			const int    navg,
			double*      csum,
			const double msum)
/* ------------------------------------------------------------------------- *
 * Subtract off mean*mean, normalise on zero-lag value.
 * ------------------------------------------------------------------------- */
{
  register int    i;
  register double mean2, norm;

  mean2 = msum * msum / navg;

  for (i = 0; i < nlag; i++) {
    csum[i] -= mean2;
  }
  norm  = csum[0];
  for (i = 0; i < nlag; i++) {
    csum[i] /= norm;
  }
}


static void printup (FILE*         fp  ,
		     const int     nlag,
		     const double* data)
/* ------------------------------------------------------------------------- *
 * Write out the autocorrelation data.
 * ------------------------------------------------------------------------- */
{
  register int i;

  for (i = 0; i < nlag; i++) fprintf (fp, "%5d %10.6g\n", i, data[i]);
}
