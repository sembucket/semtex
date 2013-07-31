/*****************************************************************************
 * asp.c: compute autospectra.
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

static void  getargs    (int, char*[], char**, int*,
			 real*, char*, int*, int*, int*);
static real  setmask    (real*, char*, int);
static int   moreinput  (FILE*);
static void  datawindow (complex*, real*, char*, int);
static void  sum        (complex*, real*, int);
static real  norm       (real*, real, real, int, int);
static void  printup    (FILE*, real*, char*, int, int, int, 
			 real,  real,  real,  int, int, int);


int main (int argc, char *argv[])
/* ------------------------------------------------------------------------- *
 * Produce autospectra from ASCII input.  Input is read in blocks of
 * blocksize points until exhausted.  If the final read is not a full one
 * it is thrown away.  If not even one full read can be done, program
 * halts.  Method of segment-averaged periodograms used.
 *
 *  Command-line arguments: These are all optional.
 *  -b <blocksize>: Segment averaging blocksize.  Default is 512.
 *      Must be an integer power of 2.  
 *  -v: Verbose output (with 8-line header).
 *  -m: Switch off mean removal, otherwise carried out on each segment
 *      before FFT.
 *  -o: Switch off segment averaging.
 *  -r <rate>:  Data sampling rate of.  Default rate is 1 Hz.
 *      <rate> is real or integer.
 *  -w <window>:    Specify a data window for leakage reduction.
 *      hann:  Von Hann window (cosine bell) Default.
 *      cost:  Cosine tip taper over 10% at each end of data.
 *      gauss: Gaussian shape.  Possible advantage for sine-wave spectra.
 *      none:  No window (Boxcar).
 *
 *  Length of records:
 *  blocksize:
 *      This is how many data points are used for each periodogram.
 *  blocksize_2:  
 *      Since the data are real, we can pack blocksize data points
 *      into blocksize_2 long complex buffer, FFT & unpack inside
 *      this buffer (since the "top" end of the buffer is redundant
 *      anyway because the data are real, there is no need to use it).
 *  blocksize_4:
 *      The default option is to do overlap segment averaging of
 *      periodograms, which results in the least variance of spectral
 *      estimate per data point recorded.  In that case, at each bite
 *      of the data file we chew off blocksize_2 points, putting
 *      them into the top half of the input data buffer, after
 *      shifting the data which were there into the lower half of the
 *      input buffer.  The input buffer is however complex, so there
 *      are blocksize_4 complex (blocksize_2 real) read in at a time.
 * ------------------------------------------------------------------------- */
{
  char     *session, window[STR_MAX];
  FILE     *fp_in, *fp_out=stdout;
  complex  *workspace, *Wtab, *inbuf;
  real     *windowmask, *autobuf;
  register  i, k;
  int       npts, navg, blocksize, blocksize_4, blocksize_2;
  real      mean, meanhat, samprate, Wss, variance;
  int       demean, overlap, verbose, starter;

  blocksize  = 512;
  samprate   = 1.0;
  demean     = 1;
  overlap    = 1;
  verbose    = 0;
  starter    = 1;
  navg       = 0;
  npts       = 0;
  meanhat    = 0.0;
  session    = 0;
  strcpy (window, "hann");

  getargs (argc, argv, &session,
	   &blocksize, &samprate, window, &demean, &overlap, &verbose);

  fp_in = (session) ? fopen (session, "r") : stdin;

  blocksize_2 = blocksize >> 1;
  blocksize_4 = blocksize >> 2;

  workspace  = cvector (0, blocksize_2 - 1);
  Wtab       = cvector (0, blocksize_2 - 1);
  inbuf      = cvector (0, blocksize_4 - 1);
  windowmask = rvector (0, blocksize   - 1);
  autobuf    = rvector (0, blocksize_2    );

  for (i = 0; i <= blocksize_2; i++) autobuf[i] = 0.0;
  Wss = setmask (windowmask, window, blocksize);
  preFFT (Wtab, blocksize_2, -1);

  if (overlap)		/* -- Do a startup half-read for overlap averaging. */
    for (i = 0; i < blocksize_4; i++)
      if (fscanf (fp_in, "%lf %lf ",
		  &inbuf[i].Re, &inbuf[i].Im) != 2)
	message ("asp", "insufficient data to half-fill buffer", ERROR);

  while (moreinput (fp_in)) {       /* -- Main processing loop. */
    if (overlap) {
      for (i = 0, mean = 0.0; i < blocksize_4; i++) {
	k = i + blocksize_4;
	workspace[i] = inbuf[i];
	if (fscanf (fp_in, "%lf %lf ",
		    &workspace[k].Re, &workspace[k].Im) != 2) break;
	mean += workspace[i].Re + workspace[i].Im +
	        workspace[k].Re + workspace[k].Im;
	npts += 4;
	inbuf[i] = workspace[k];
      }
      if ((i != blocksize_4) && starter)
	message ("asp", "insufficient data to fill buffer", ERROR);
    } else {			/* -- No overlap. */
      for (i = 0, mean = 0.0; i < blocksize_2; i++) {
	if (fscanf (fp_in, "%lf %lf ",
		    &workspace[i].Re, &workspace[i].Im) != 2) break;
	mean += workspace[i].Re + workspace[i].Im;
	npts += 2;
      }
      if ((i != blocksize_2) && starter)
	message ("asp", "insufficient data to fill buffer", ERROR);
    }

    starter = 0;
    mean    = mean / blocksize;
    meanhat = meanhat + mean;
    if (demean)
      for (i = 0; i < blocksize_2; i++) {
	workspace[i].Re -= mean;
	workspace[i].Im -= mean;
      }
    if (!strstr (window, "none"))
      datawindow (workspace, windowmask, window, blocksize_2);
    rcFFT (workspace, blocksize_2, Wtab, blocksize_2, 1);
    sum   (workspace, autobuf, blocksize_2);
    navg++;
  }

  meanhat  = meanhat / navg;
  variance = norm (autobuf, Wss, samprate, blocksize_2, navg);

  printup (fp_out, autobuf, window, blocksize, navg, npts,
	   samprate, meanhat, variance, demean, overlap, verbose);

  return EXIT_SUCCESS;
}


static void getargs (int    argc      ,
		     char*  argv[]    ,
		     char** session   ,
		     int*   blocksize ,
		     real*  samprate  ,
		     char*  window    ,
		     int*   removemean,
		     int*   overlap   ,
		     int*   verbose   )
/* ------------------------------------------------------------------------- *
 * Process command-line arguments.
 * ------------------------------------------------------------------------- */
{
  char   c;
  static char *usage = 
    "Usage: asp [options] [file]\n"
    "  options are:\n"
    "  -h             ... print this message\n"
    "  -b blocksize   ... size of data segments         [D: 512]\n"
    "  -m             ... switch off mean removal\n"
    "  -o             ... switch off overlap averaging\n"
    "  -r rate        ... sampling rate                 [D: 1Hz]\n"
    "  -v             ... verbose output\n"
    "  -w window      ... one of \"none\", \"hann\", \"gauss\", \"cost\"\n";
  
  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fprintf(stderr, usage);
      exit(0);
      break;
    case 'b':
      if (*++argv[0])
	*blocksize = atoi(*argv);
      else {
	--argc;
	*blocksize = atoi(*++argv);
      }
      if (!ispow2(*blocksize))
	message("asp", " blocksize must be a power of 2", ERROR);
      break;
    case 'm':
      *removemean = 0;
      break;
    case 'o':
      *overlap = 0;
      break;
    case 'r':
      if (*++argv[0])
	*samprate = atof(*argv);
      else {
	--argc;
	*samprate = atof(*++argv);
      }
      break;
    case 'v':
      *verbose = 1;
      break;
    case 'w':
      if(*++argv[0])
	strcpy(window, *argv);
      else {
	--argc;
	strcpy(window, *++argv);
      }
      if (!(   strstr(window, "none")
	    || strstr(window, "hann")
	    || strstr(window, "gauss")
	    || strstr(window, "cost")))
	message ("asp: unknown window specifier", window, ERROR);
      break;
    default:
      fprintf(stderr, "asp: illegal option: %c\n", c);
      break;
    }
  
  if (argc == 1) *session = *argv;
}


static real setmask (real* mask  ,
		     char* window,
		     int   block )
/* ------------------------------------------------------------------------- *
 * Set up a mask for data-windowing.
 * ------------------------------------------------------------------------- */
{
  int  j, block_2, tenpercent;
  real twopi_block, Wss;

  Wss     = 0.0;
  block_2 = block >> 1;

  if (strstr(window, "hann")) {
    twopi_block = (M_PI + M_PI) / block;
    for (j = 0; j < block; j++) {
      mask[j] = 0.5 * (1.0 - cos(twopi_block*j));
      Wss += SQR(mask[j]);
    }
    Wss *= block;
  } else if (strstr(window, "cost")) {
    tenpercent = block / 10;
    for (j = 0; j < tenpercent; j++) {
      mask[j]  = 0.5 * (1.0 - cos(j * M_PI/tenpercent));
      Wss     += 2.0 * SQR(mask[j]);
    }
    Wss = block * (Wss + block - 2 * tenpercent);
  } else if (strstr(window, "gauss")) {
    for (j = 0; j < block; j++) {
      mask[j] = exp( -0.5 * SQR(6.0*(j - block_2)/block) );
      Wss += SQR(mask[j]);
    }
    Wss *= block;
  } else /* -- No window, uniform weighting, masking non carried out. */
    Wss = block;

  return Wss;
}


static int moreinput (FILE *fp)
/* ------------------------------------------------------------------------- *
 * Check ASCII input file.
 * ------------------------------------------------------------------------- */
{
  register c;

  while (isspace (c = fgetc (fp))) ;
  
  if (c == EOF) {          return 0;
  } else { ungetc (c, fp); return 1;
  }
}


static void datawindow (complex*  work   ,
			real*     mask   ,
			char*     window ,
			int       nyquist)
/* ------------------------------------------------------------------------- *
 * Apply data window for leakage reduction.
 * ------------------------------------------------------------------------- */
{
  register i, tenpercent;

  if (strstr (window, "cost")) {
    tenpercent = 2 * nyquist / 10;
    work[0].Re *= mask[0];
    for (i = 1; i < tenpercent; i++) {
      if (i % 2) {
	work[          (i-1)/2].Im *= mask[i];
	work[nyquist - (i+1)/2].Im *= mask[i];
      } else {
	work[              i/2].Re *= mask[i];
	work[nyquist     - i/2].Re *= mask[i];
      }
    }
  } else {
    for (i = 0; i < nyquist; i++) {
      work[i].Re *= mask[i + i];
      work[i].Im *= mask[i + i + 1];
    }
  }
}


static void sum (complex* work   ,
		 real*    autobuf,
		 int      nyquist)
/* ------------------------------------------------------------------------- *
 * Add contribution to averaging buffer.
 * ------------------------------------------------------------------------- */
{
  register int i;

  autobuf[0]       +=       SQR(work[0].Re);
  autobuf[nyquist] += 2.0 * SQR(work[0].Im);
  for (i = 1; i < nyquist; i++)
    autobuf[i] += 2.0 * (SQR(work[i].Re) + SQR(work[i].Im));
}



static real norm (real* autobuf,
		  real  Wss    , 
		  real  rate   ,
		  int   nyquist,
		  int   navg   )
/* ------------------------------------------------------------------------- *
 * Normalize spectrum, return variance.  Power scaling (see Press et al.).
 * ------------------------------------------------------------------------- */
{
  real     variance = 0.0, scale, deltaF;
  register i;

  deltaF   = 0.5 *  rate / nyquist;
  scale    = 1.0 / (Wss * navg * deltaF);

  for (i = 0; i <= nyquist; i++) {
    autobuf[i] *= scale;
    variance   += autobuf[i];
  }
  
  return deltaF * variance;
}


static void printup (FILE* fp       ,
		     real* autobuf  ,
		     char* window   ,
		     int   blocksize, 
		     int   navg     ,
		     int   npts     ,
		     real  samprate ,
		     real  meanhat  ,
		     real  variance ,
		     int   demean   ,
		     int   overlap  ,
		     int   verbose  )
/* ------------------------------------------------------------------------- *
 * Print up results.
 * ------------------------------------------------------------------------- */
{
  static char *hdr_fmt[] = {
    "%-25s "      "%s\n",
    "%-25s "      "Created\n",
    "%-25d "      "Number of data\n",
    "%-25d "      "Transform length\n",
    "%-25d "      "Ensemble averages\n",
    "%-25s "      "Data window\n",
    "%-25.6g "    "Sampling rate\n",
    "%-25d "      "Frequency bins\n",
    "%-25.6g "    "Mean value\n",
    "%-25.6g "    "Variance\n"
    };
  static  char name[] = "asp";
  time_t  tp;
  int     i, bs_2 = blocksize >> 1;
  real    deltaF = 0.5 * samprate / bs_2, F = 0.0;
  char    buf[BUFSIZ];

  if (verbose) {
    sprintf (buf, "Autospectrum");
    if (demean)
      strcat (buf, " [mean removed]");
    if (overlap)
      strcat (buf, " [overlap averaged]");
    fprintf (fp, hdr_fmt[0], name, buf);
    
    tp = time ((time_t*) NULL);
    strftime (buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));
    fprintf (fp, hdr_fmt[1], buf);
    fprintf (fp, hdr_fmt[2], npts);
    fprintf (fp, hdr_fmt[3], bs_2);
    fprintf (fp, hdr_fmt[4], navg);
    fprintf (fp, hdr_fmt[5], window);
    fprintf (fp, hdr_fmt[6], samprate);
    fprintf (fp, hdr_fmt[7], bs_2 + 1);
    fprintf (fp, hdr_fmt[8], meanhat);
    fprintf (fp, hdr_fmt[9], variance);
    fprintf (fp, "\n");
  }

  for (i = 0; i <= bs_2; i++) {
    fprintf (fp, "%-16g %-16g\n", F, autobuf[i]);
    F += deltaF;
  }
}
