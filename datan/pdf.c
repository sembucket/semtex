/*****************************************************************************
 * PDF: produce probability density function from columnar format input.
 *
 * Synopsis:
 * ---------
 * Read ASCII data from stdin or named file, place in temporary scratch file
 * and record minumum & maximum values, number of data.  Create an estimate
 * of the pdf by accumulating sums of relative frequencies.  Optional
 * command-line argument supplies number of bins in which sums are accumul-
 * ated.  Print pdf estimates to stdout or optional named file, together
 * with an optional header (specified using the command-line argument -v)
 * which contains information including moments (mean, variance...) of the
 * data, up to fourth order (kurtosis).
 *
 * Usage:
 * ------
 * pdf [-h] [-n <#bins>] [-o <output>] [-v] [input]
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "FFTutil.h"

static void getargs   (int, char*[], int*, int*, FILE**, FILE**);
static void refile    (FILE*, FILE*, int*, real*, real*, real*);
static void accum     (FILE*, real*, real, real, int, int, real,
		       real*, real*, real*, real*, real*);
static void printhead (FILE*, real,   real,  real,  real, int, int,
                       real,  real, real, real,  real, real);
static void printup   (FILE*, real*, real, real, int);

#define    DEFBIN     21
#define    TINY       1.0E-30


int main(int argc, char *argv[])
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  real *pdf;
  FILE  *fp_Store,
        *fp_in   = stdin,
        *fp_out  = stdout;
  int    Verbose = 0;

  int    i,
         Npts  = 0,
         NBins = DEFBIN;

  real  BinWidth,
         Min      = 0.0,
         Max      = 0.0,
         PeakVal  = 0.0,
         Mean     = 0.0,
         SDev     = 0.0,
         AbsDev   = 0.0,
         Variance = 0.0,
         Skewness = 0.0,
         Kurtosis = 0.0;


  getargs (argc, argv, &NBins, &Verbose, &fp_in, &fp_out);

  pdf = rvector (0, NBins-1);

  fp_Store = tmpfile ();
  refile (fp_in, fp_Store, &Npts, &Min, &Max, &Mean);
  BinWidth = (Max - Min) / NBins;
  if (BinWidth < TINY) message ("pdf", "bin width too small", ERROR);

  rewind (fp_Store);
  accum  (fp_Store, pdf, Min, BinWidth, Npts, NBins,
	  Mean, &AbsDev, &SDev, &Variance, &Skewness, &Kurtosis);

  for (i = 0; i < NBins; i++) {
    pdf[i] /= BinWidth * Npts;
    PeakVal = MAX(pdf[i], PeakVal);
  }

  if (Verbose)
    printhead (fp_out,  Min, Max, PeakVal, BinWidth, NBins, Npts,
	       Mean, AbsDev, SDev, Variance, Skewness, Kurtosis);

  printup (fp_out, pdf, BinWidth, Min, NBins);

  return EXIT_SUCCESS;
}


static void getargs (int    argc,    char  *argv[],
		     int   *NBins,   int   *Verbose,
		     FILE **fp_in,   FILE **fp_out)
/* ------------------------------------------------------------------------- *
 * Parse command line.
 * ------------------------------------------------------------------------- */
{
  char c;
  char line[FILENAME_MAX];
  static char usage[] = 
    "usage: pdf [options] [input]\n"
    "[options]:\n"
    "-h          ... display this message\n"      
    "-n <#bins>  ... create PDF with #bins bins          [Default: 21]\n"
    "-o <output> ... write output to named file          [Default: stdout]\n"
    "-v          ... set verbose (add header to output)\n";


  while (--argc && (*++argv)[0] == '-') {
    switch (c = *++argv[0]) {
    case 'h':
      fputs(usage, stderr);
      exit(0);
      break;
    case 'n':
       if (*++argv[0])
	*NBins = atoi(*argv);
      else {
	*NBins = atoi(*++argv);
	argc--;
      }
      break;
    case 'o':
      if (*++argv[0])
	*fp_out = fopen(*argv, "w");
      else {
	*fp_out = fopen(*++argv, "w");
	argc--;
      }
      break;
    case 'v':
      *Verbose = 1;
      break;

    default:
      message ("pdf", "unknown option", WARNING);
      fprintf (stderr, usage);
      exit (EXIT_FAILURE);
      break;
    }
  }

  if (argc == 1)
    if ((*fp_in = fopen (*argv, "r")) == (FILE*) NULL)
      message ("pdf", "couldn't open input file", ERROR);

}


static void refile (FILE *fp_in, FILE *fp_tmp,
		    int *npts, real *Min, real *Max, real *Mean)
/* ------------------------------------------------------------------------- *
 * Fill temp file with data, getting parameters needed for creation of PDF.
 * ------------------------------------------------------------------------- */
{
  real      datum;
  char      line[FILENAME_MAX];
  const int DP = sizeof (real) == sizeof (double);

  while (fgets (line, FILENAME_MAX, fp_in)) {
    if   (DP) sscanf (line, "%lf", &datum);
    else      sscanf (line, "%f",  &datum);
    if (!*npts)
      *Min = (*Max = datum);
    else {
      *Min = MIN (datum, *Min);
      *Max = MAX (datum, *Max);
    }
    *Mean += datum;
    ++*npts;
    fwrite (&datum, sizeof(real), 1, fp_tmp);
  }

  *Mean /= *npts;
}





static void accum(FILE  *fp,
		  real *pdf,
		  real  min,  real binwidth,
		  int    npts, int nbins, 
		  real  mean,
		  real *AbsDev,
		  real *Sdev,
		  real *Vari,
		  real *Skew,
		  real *Kurt)
/*===========================================================================*
 * Make all the moments.                                                     *
 *===========================================================================*/
{
  real datum, dev, prod;
  int   bin;


  while (!feof(fp)) {
    fread(&datum, sizeof(real), 1, fp);
    dev     = datum - mean;
    *AbsDev = *AbsDev + fabs(dev);
    prod    =  dev * dev;
    *Vari   = *Vari + prod;
    prod    =  prod * dev;
    *Skew   = *Skew + prod;
    prod    =  prod * dev;
    *Kurt   = *Kurt + prod;
    bin     = ((datum - min) / binwidth);
    if (bin > nbins-1) bin = nbins - 1;
    else if (bin < 0)  bin = 0;
    pdf[bin] += 1.0;
  }
    
  *AbsDev = *AbsDev /  npts;
  *Vari   = *Vari   / (npts - 1);
  *Sdev   = sqrt(*Vari);
  if (*Vari > 0.0) {
    *Skew = *Skew / (npts * SQR(*Vari) / *Sdev);
    *Kurt = *Kurt / (npts * SQR(*Vari)) - 3.0;
  }
}





static void printhead(FILE  *fp,
		      real  min,   real max,  real pk,  real bnwd, 
		      int    nbins, int   npts,
		      real  mean,  real adev, real sdev,
		      real  vari,  real skew, real kurt )
/*===========================================================================*
 * Write header info.                                                        *
 *===========================================================================*/
{
  static char *hdr_fmt[] = {
    "%-25s "      "Probability density function\n",
    "%-25s "      "Created\n",
    "%-25d "      "Number of data\n",
    "%-25d "      "Number of bins\n",
    "%-25.6g "    "Bin width\n",
    "%-25.6g "    "Maximum probability density\n",
    "%-25.6g "    "Minimum value\n",
    "%-25.6g "    "Maximum value\n",
    "%-25.6g "    "Mean value\n",
    "%-25.6g "    "Standard deviation\n",
    "%-25.6g "    "Absolute deviation\n",
    "%-25.6g "    "Variance\n",
    "%-25.6g "    "Skewness\n",
    "%-25.6g "    "Kurtosis\n\n"
    };
  static char name[] = "pdf";
  time_t tp;
  char   buf[BUFSIZ];


  tp = time((time_t*) NULL);
  strftime(buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));

  fprintf(fp, hdr_fmt[0],  name);
  fprintf(fp, hdr_fmt[1],  buf);
  fprintf(fp, hdr_fmt[2],  npts);
  fprintf(fp, hdr_fmt[3],  nbins);
  fprintf(fp, hdr_fmt[4],  bnwd);
  fprintf(fp, hdr_fmt[5],  pk);
  fprintf(fp, hdr_fmt[6],  min);
  fprintf(fp, hdr_fmt[7],  max);
  fprintf(fp, hdr_fmt[8],  mean);
  fprintf(fp, hdr_fmt[9],  sdev);
  fprintf(fp, hdr_fmt[10], adev);
  fprintf(fp, hdr_fmt[11], vari);
  fprintf(fp, hdr_fmt[12], skew);
  fprintf(fp, hdr_fmt[13], kurt);
}





static void printup(FILE  *fp,
		    real *pdf,
		    real  bnwd,
		    real  min,
		    int    nbins )
/*===========================================================================*
 * Write out pdf data.                                                       *
 *===========================================================================*/
{
  register i;


  for (i=0; i<nbins; i++)
    fprintf(fp, "%10.6g %10.6g\n", min+(i+0.5)*bnwd, *(pdf + i));
}
