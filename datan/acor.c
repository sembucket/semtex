/*****************************************************************************
 * ACOR: produce autocorrelation function from columnar format input.
 *                                                                           *
 * Synopsis:
 * ---------
 * Read ASCII data from stdin or named file, place in temporary scratch file
 * and record mean value, number of data.  Optionally remove the mean.
 * Zero-pad by the number of lags specified on the command line (to avoid
 * overlap effects) and then, if required, pad with more zeros to produce a
 * number of data points which is an integer power of two (for FFT).
 * Compute autocorrelation using FFT methods.  Optionally nondimensionalize
 * by dividing through by mean-squared value.  Print out to stdout (or
 * optionally a named file), with or without an optional header.
 *
 * Usage:
 * ------
 * acor [-d] [-h] [-m] [-n <#lags>] [-o <output>] [-r <rate>] [-v] [input]
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "FFTutil.h"

static void getargs  (int, char*[],int*,int*,int*,int*,real*,FILE**,FILE**);
static void refile   (FILE*, FILE*, int*, real*);
static void pad      (real*, FILE*, int, real, int, int);
static void nondimen (real*, int);
static void normalize(real*, int, int);
static void printhead(FILE*, int, int, int, int, int, int, real, int, real);
static void printup  (FILE*, const real*, int, real);

#define DEFLAG   32
#define FORWARD  1
#define INVERSE  0


int main (int argc, char *argv[])
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  complex      *Zbuf, *Wtab;
  FILE         *fp_Store,
               *fp_in  = stdin,
               *fp_out = stdout;
  int           Verbose = 0, NonDim = 0, DeMean = 0, NData = 0, NLag = DEFLAG;
  int           NPad, TranLen;
  real          SampRate = 1.0, Mean = 0.0;
  register int  i;

  /* -- Parse command line, place data in temporary file. */

  getargs (argc,argv,&NLag,&NonDim,&Verbose,&DeMean,&SampRate,&fp_in,&fp_out);

  fp_Store = tmpfile();
  refile (fp_in, fp_Store, &NData, &Mean);

  /* -- Set up storage, precompute FFT angular factors. */

  TranLen = roundpow2 (NData + NLag) >> 1;
  NPad    = (TranLen << 1) - NData;
  
  Zbuf = cvector(0, TranLen-1);
  Wtab = cvector(0, TranLen-1);
  
  preFFT (Wtab, TranLen, -1);
    
  /* -- Do padding to make up power of 2 if needed & to avoid overlap. */

  rewind (fp_Store);
  pad    ((real*)Zbuf, fp_Store, DeMean, Mean, NData, NPad);

  /* -- Make autocorrelation by backtransform of product of DFT & conjugate. */
  
  rcFFT (Zbuf, TranLen, Wtab, TranLen, FORWARD);
  
  Zbuf[0].Re = SQR (Zbuf[0].Re);
  Zbuf[0].Im = SQR (Zbuf[0].Im);
  
  for (i = 1; i < TranLen; i++) {
    Zbuf[i].Re = SQR (Zbuf[i].Re) + SQR (Zbuf[i].Im);
    Zbuf[i].Im = 0.0;
  }
   
  rcFFT (Zbuf, TranLen, Wtab, TranLen, INVERSE);
  
  /* -- Scale, print, exit. */

  if (NonDim)
    nondimen ((real*)Zbuf, TranLen << 1);
  else
    normalize ((real*)Zbuf, TranLen, NData);
    
  if (Verbose)
    printhead (fp_out,
	       Verbose, NData, NPad, NLag, TranLen, DeMean, Mean, NonDim,
	       SampRate);
 
  printup (fp_out, (real*)Zbuf, NLag, SampRate);

  return EXIT_SUCCESS;
}


static void getargs (int    argc,    char**  argv,
		     int*   NLag,    int*    NonDim,
		     int*   Verbose, int*    DeMean,
		     real*  SampRate,
		     FILE** fp_in,   FILE**  fp_out)
/* ------------------------------------------------------------------------- *
 * Parse command line.
 * ------------------------------------------------------------------------- */
{
  char c;
  char line[FILENAME_MAX];
  static char usage[] = 
    "usage: acor [options] [input]\n"
    "[options]:\n"
    "-d          ... create dimensionless autocorrelation\n"
    "-h          ... display this message\n"
    "-m          ... remove mean value before computation\n"
    "-n <#lags>  ... compute this many lags               [Default: 32]\n"
    "-o <output> ... write output to named file           [Default: stdout]\n"
    "-r <rate>   ... Sampling rate [Hz]                   [Default: 1 Hz]\n"
    "-v          ... set verbose (add header to output)\n";


  while (--argc && (*++argv)[0] == '-') {
    switch (c = *++argv[0]) {
    case 'd':
      *NonDim = 1;
      break;
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;
    case 'm':
      *DeMean = 1;
      break;
    case 'n':
      if (*++argv[0])
	*NLag = atoi (*argv);
      else {
        *NLag = atoi (*++argv);
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
    case 'r':
      if (*++argv[0])
    	*SampRate = atof (*argv);
      else {
    	*SampRate = atof (*++argv);
    	argc--;
      }
      break;
    case 'v':
      *Verbose = 1;
      break;

    default:
      message ("acor", "unknown option", WARNING);
      break;
    }
  }

  if (argc == 1)
    if ((*fp_in = fopen (*argv, "r")) == (FILE*) NULL)
      message ("acor", "couldn't open input file", ERROR);

}


static void refile (FILE*  fp_in, FILE*  fp_tmp,
		    int*   npts,  real*  Mean  )
/* ------------------------------------------------------------------------- *
 * Fill temp file with data.
 * ------------------------------------------------------------------------- */
{
  real  datum;
  char  line[FILENAME_MAX];

  while (fgets (line, FILENAME_MAX, fp_in)) {  	
    sscanf (line, "%f", &datum);
    *Mean += datum;
    ++*npts;
    fwrite (&datum, sizeof(real), 1, fp_tmp);
  }

  *Mean /= *npts;
}


static void pad (real*  Data,   FILE* fp,
		 int    DeMean, real  Mean, int NData, int NPad)
/* ------------------------------------------------------------------------- *
 * Recall data from temporary file, remove mean if applicable, zero pad.
 * ------------------------------------------------------------------------- */
{
  register int i;
  
  for (i = 0; i < NData; i++) {
    fread(Data + i, sizeof(real), 1, fp);
    if (DeMean) *(Data + i) -= Mean;
  }
  
  for (i = NData; i < NData+NPad; i++) *(Data + i) = 0.0;
}


static void nondimen (real*  Data,  int  N)
/* ------------------------------------------------------------------------- *
 * Scale to make first (zero lag) value = 1.
 * ------------------------------------------------------------------------- */
{
  register int    i;
  register real   Var = *Data;
    
  for (i = 0; i < N; i++)  *(Data + i) /= Var;
}


static void normalize (real* Data, int N, int NData)
/* ------------------------------------------------------------------------- *
 * Scale to make first (zero lag) value = mean squared value.
 * ------------------------------------------------------------------------- */
{
  register int  i;
  const    real    Fac = 1.0 / (N * NData);

  for (i = 0; i < N; i++) *(Data + i) *= Fac;
}


static void printhead (FILE  *fp      ,
		       int    Verbose ,
		       int    NData   ,
		       int    NPad    ,
		       int    NLag    ,
		       int    TranLen ,
		       int    DeMean  ,
		       real   Mean    ,
		       int    NonDim  ,
		       real   SampRate)
/* ------------------------------------------------------------------------- *
 * Write header information.
 * ------------------------------------------------------------------------- */
{
  static char *hdr_fmt[] = {
    "%-25s "      "%s\n",
    "%-25s "      "Created\n",
    "%-25d "      "Number of data\n",
    "%-25d "      "Number of lags\n",
    "%-25d "      "Zero padding\n",
    "%-25d "      "Transform length\n",
    "%-25.6g "    "%s\n",
    "%-25.6g "    "Sampling rate\n"
    };
  static char name[] = "acor";
  time_t tp;
  char   buf[BUFSIZ];

  sprintf (buf, "Autocorrelation");
  if (NonDim)
    strcat (buf, " [dimensionless]");
  fprintf (fp, hdr_fmt[0], name, buf);

  tp = time ((time_t*) NULL);
  strftime (buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));
  fprintf (fp, hdr_fmt[1], buf);

  fprintf (fp, hdr_fmt[2], NData);
  fprintf (fp, hdr_fmt[3], NLag);
  fprintf (fp, hdr_fmt[4], NPad);
  fprintf (fp, hdr_fmt[5], TranLen);
	   
  sprintf (buf, "Mean value");
  if (DeMean)
    strcat (buf, " [removed]");
  fprintf (fp, hdr_fmt[6], Mean, buf);
  fprintf (fp, hdr_fmt[7], SampRate);
  fprintf (fp, "\n");
}


static void printup (FILE* fp, const real* Data, int NLag, real SampRate)
/* ------------------------------------------------------------------------- *
 * Write out the autocorrelation data.
 * ------------------------------------------------------------------------- */
{
  register int i;

  for (i = 0; i < NLag; i++)
    fprintf (fp, "%10.6g %10.6g\n", i/SampRate, Data[i]);
}
