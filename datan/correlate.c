/*****************************************************************************
 * correlate.c: compute cross-correlation of two columns of input data.
 *
 * Synopsis
 * --------
 * Produce cross-correlation from 2 columns of ASCII input.  Input is
 * read in chuncks of blocksize points until exhausted.  If the final
 * read is not a full one it is thrown away.  If not even one full
 * read can be done, program halts.
 *
 * Usage
 * -----
 * correlate [options] [file]
 *  -b <blocksize>: Segment averaging blocksize.  Default is 512, must be even.
 *  -k <lag>: Number of lags (+/-) over which to compute correlation.
 *  -v: Verbose output (with 8-line header).
 *  -o: Switch off segment averaging.
 *  -r <rate>:  Data sampling rate.  Default rate is 1 (Hz).
 *
 * Use routine correl_() from Numerical Recipes 2e, but compute the
 * "true" (dimensionless) cross correlation
 *
 *   R[u(t),v(t+lag)] = E[(u-E(u))*(v-E(v))]/[sdev(u)*sdev(v)]
 *
 * NBNBNBNB: Owing to a limitation in correl_() the maximum FFT size
 * allowed is 2^16 i.e. 65536 (as originally published this value was
 * 4096).
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include "FFTutil.h"

static char prog[]  = "correlate";
static int  DP      = 0;	/* Global precision flag. */
static int  verbose = 0;

#if 1
#define NVREG  8
#define BLKMAX 65536

int    _alpIreg[NVREG];		/* For FORTRAN linkage. */
char   _alpCreg[NVREG];
float  _alpSreg[NVREG];
double _alpDreg[NVREG];

void    correl_ (real*, real*, int*, real*);
#define correl(a, b, n, c) (_alpIreg[0]=n, correl_(a, b, _alpIreg, c))
#endif

static void  getargs   (int, char**, char**, int*, int*, real*, int*);
static int   moreinput (FILE*);
static void  preproc   (real*, real*, int);
static void  printup   (FILE*, char*, real*, int, int, int, real, int);


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  char         *session = 0;
  FILE         *fp_in, *fp_out=stdout;
  real         *work, *buffa, *buffb, *buffc, *average;
  int          nlag, npts, navg, blocksize, blocksize_2, bufflen, npad;
  real         samprate, norm;
  int          overlap, starter;
  register int i, k;
  
  DP = (sizeof (real) == sizeof (double) ) ? 1 : 0;

  blocksize  = 512;
  nlag       = 255;
  samprate   = 1.0;
  overlap    = 1;
  starter    = 1;
  navg       = 0;

  getargs (argc, argv, &session, &blocksize, &nlag, &samprate, &overlap);

  fp_in = (session) ? fopen (session, "r") : stdin;

  blocksize_2 = blocksize >> 1;
  bufflen     = roundpow2 (blocksize + nlag);
  npad        = bufflen - blocksize; 

  if (bufflen > BLKMAX)
    message (prog, "exceeded the maximum allowed buffer length", ERROR);

  work        = rvector (0, 6 * bufflen - 1);
  buffa       = work;
  buffb       = buffa   + bufflen;
  average     = buffb   + bufflen;
  buffc       = average + bufflen;

  for (i = 0; i < bufflen; i++) average[i] = 0.0;

  if (overlap) {
    for (i = blocksize_2; i < blocksize; i++)
      if (fscanf (fp_in, (DP) ? "%lf %lf" : "%f %f", buffa + i, buffb + i)
	  != 2) break;
    if (i != blocksize)
      message (prog, " insufficient data", ERROR);
  }

  while (moreinput (fp_in)) {       /* -- Main processing loop. */
    if (overlap)
      for (npts = 0, i = blocksize_2; i < blocksize; i++, npts += 2) {
	k = i - blocksize_2;
	buffa[k] = buffa[i];
	buffb[k] = buffb[i];
	if (fscanf (fp_in, (DP) ? "%lf %lf" : "%f %f", buffa + i, buffb + i)
	    != 2) break;
      }
    else
      for (npts = 0, i = 0; i < blocksize; i++, npts++)
	if (fscanf (fp_in, (DP) ? "%lf %lf": "%f %f", buffa + i, buffb + i)
	    != 2) break;

    if (npts != blocksize)
      if (starter) message (prog, " insufficient data", ERROR);
      else break;
    if (starter) starter = 0;

    for (i = 0; i < npad; i++) buffa[blocksize+i] = buffb[blocksize+i] = 0.0;

    preproc (buffa, buffb, blocksize);

    correl (buffa, buffb, bufflen, buffc);

    for (i = 0; i < bufflen; i++) average[i] += buffc[i];
    navg++;
  }

  norm = 1.0 / (navg * (blocksize - 1));
  for (i = 0; i < bufflen; i++) average[i] *= norm;

  printup (fp_out, session, average, bufflen, nlag, navg, samprate, overlap);

  return EXIT_SUCCESS;
}


static void getargs (int    argc     ,
		     char** argv     ,
		     char** session  ,
		     int*   blocksize,
		     int*   nlag     ,
		     real*  samprate ,
		     int*   overlap  )
/* ------------------------------------------------------------------------- *
 * Process command-line arguments.
 * ------------------------------------------------------------------------- */
{
  char c;
  static char *usage = 
    "Usage: correlate [options] [file]\n"
    "options:\n"
    "-h           ... print this message\n"
    "-b blocksize ... size of data segments         [D: 512]\n"
    "-k lags      ... number of lags                [D: 255]\n"
    "-o           ... switch off overlap averaging\n"
    "-r rate      ... sampling rate                 [D: 1Hz]\n"
    "-v           ... verbose output\n";
  
  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fprintf(stderr, usage);
      exit(0);
      break;
    case 'b':
      if (*++argv[0]) *blocksize = atoi(*argv);
      else {--argc; *blocksize = atoi(*++argv);}
      if (*blocksize & 1) message (prog, " blocksize must be even", ERROR);
      break;
    case 'k':
      if (*++argv[0]) *nlag = atoi(*argv);
      else {--argc; *nlag = atoi(*++argv);}
      break;
    case 'o':
      *overlap = 0;
      break;
    case 'r':
      if (*++argv[0]) *samprate = atof(*argv);
      else {--argc; *samprate = atof(*++argv);}
      break;
    case 'v':
      verbose = 1;
      break;
    default:
      fprintf (stderr, "correlate: illegal option: %c\n", c);
      break;
    }

  if (*blocksize <= *nlag + *nlag + 1)
    message (prog, " number of lags must be < 0.5 blocksize", ERROR);
  if (argc == 1) *session = *argv;
}


static int moreinput (FILE *fp)
/* ------------------------------------------------------------------------- *
 * Check ASCII input file.
 * ------------------------------------------------------------------------- */
{
  register int c;
  while (isspace(c = getc(fp)));
  if (c == EOF) {return 0;} else {ungetc(c, fp); return 1;}
}


static void printup (FILE* fp      ,
		     char* session ,
		     real* average ,
		     int   bufflen , 
		     int   nlag    ,
		     int   navg    ,
		     real  samprate,
		     int   overlap )
/* ------------------------------------------------------------------------- *
 * Print up results.
 * ------------------------------------------------------------------------- */
{
  static char *hdr_fmt[] = {
    "%-25s "      "%s\n",
    "%-25s "      "Created\n",
    "%-25d "      "Number of lags\n",
    "%-25d "      "Ensemble averages\n",
    "%-25.6g "    "Sampling rate\n"
    };
  time_t tp;
  int    i, k;
  real   deltaT = 1.0 / samprate, T;
  char   buf[STR_MAX];

  if (verbose) {
    sprintf (buf, "Cross correlation");
    if (overlap) strcat (buf, " [overlap averaged]");
    fprintf (fp, hdr_fmt[0], session, buf);
    
    tp = time ((time_t*) NULL);
    strftime (buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));
    fprintf (fp, hdr_fmt[1], buf);
    fprintf (fp, hdr_fmt[2], nlag);
    fprintf (fp, hdr_fmt[3], navg);
    fprintf (fp, hdr_fmt[4], samprate);
  }

  for (i = nlag; i; i--)
    fprintf (fp, "%-16g %-16g\n", i * -deltaT, average[bufflen - i]);

  for (i = 0; i <= nlag; i++)
    fprintf (fp, "%-16g %-16g\n", i *  deltaT, average[i]);
}


static void  preproc (real* a,
		      real* b,
		      int   N)
/* ------------------------------------------------------------------------- *
 * Subtract mean from each process and normalise by its standard deviation.
 * ------------------------------------------------------------------------- */
{
  int  i;
  real varia = 0.0, varib = 0.0, meana = 0.0, meanb = 0.0;

  for (i = 0; i < N; i++) {
    meana = meana + a[i];
    meanb = meanb + b[i];
  }
  meana = meana / N;
  meanb = meanb / N;

  for (i = 0; i < N; i++) {
    a[i]  -= meana;
    b[i]  -= meanb;
    varia += a[i] * a[i];
    varib += b[i] * b[i];
  }
  varia = varia / (N-1);	/* -- Unbiased estimates of variance. */
  varib = varib / (N-1);

  varia = 1.0 / sqrt(varia);    /* -- (Inverse) standard deviations. */
  varib = 1.0 / sqrt(varib);

  for (i = 0; i < N; i++) {
    a[i] *= varia;
    b[i] *= varib;
  }
}


