/*****************************************************************************
 *                  MEMSPEC: ALL-POLE SPECTRAL ANALYSIS.                     *
 *                                                                           *
 * Synopsis:                                                                 *
 * ---------                                                                 *
 * Given a vector of real numbers on stdin, compute spectral estimates for   *
 * the data using the Maximum Entropy Method.  Write the estimates, with     *
 * header information, to stdout.  Frequencies may be interpreted in three   *
 * ways; by default, the sampling frequency is taken to be unity, and the    *
 * estimates are scaled to fit within the Nyquist range 0..0.5 (with densit- *
 * ies scaled such that the integral from 0 to 0.5 equals mean-square value).*
 * If a sampling frequency is given on the command line, the densities are   *
 * normalized so that the integral from zero to the Nyquist frequency is the *
 * mean square value.  Finally, the frequencies may be converted to dimen-   *
 * sionless form, on the basis of velocity scale and/or length scale given   *
 * on the command line.  In this case, the scaling is done so that the area  *
 * of the spectrum from 0 to the reduced Nyquist frequency is the mean       *
 * square value.  (fr = reduced frequency = fD/U.)  The frequency range      *
 * over which the spectral estimates are computed is always interpreted in   *
 * terms of reduced frequency.                                               *
 *                                                                           *
 * Usage:                                                                    *
 * ------                                                                    *
 * memspec [-m poles -n evals -r fsamp -U vel -D dia -s start -e end] [file] *
 *                                                                           *
 * Defaults:                                                                 *
 * ---------                                                                 *
 * Number of poles:            DEFPOL;                                       *
 * Number of evaluations:      DEFEVL;                                       *
 * Maximum number of data:     MAXPNT;                                       *
 * Sampling Frequency:         1.0;                                          *
 * Velocity Scale:             1.0;                                          *
 * Length Scale:               1.0;                                          *
 *                                                                           *
 * Reference:                                                                *
 * ----------                                                                *
 * Press, W.H. et al. 1992, "Numerical Recipes", 2e, C.U.P., Section 13.7.   *
 *****************************************************************************/

/*------------------*
 * RCS Information: *
 *------------------*/
static char
  RCSid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include <nrutil.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#define  MAXPNT  8192
#define  DEFPOL  10
#define  DEFEVL  512


static void  memcof   (float data[], int n, int m, float *xms, float d[]);
static float evlmem   (float fdt, float d[], int m, float xms);
static void  getargs  (int argc, char *argv[], FILE **fp_in, 
	               float *Fstart, float *Fend, float *Fsamp,
		       float *Diam,   float *Vel,  int *Neval, int *Npoles);
static void  getdata  (FILE *fp, float *Data, int *Ndata, int MAX);
static void  printhead(FILE *fp, float Fs, float Fe, int Nd, int Np);
static void  doscale  (float  fstart, float *fend, float fsamp,
                       float  dia,    float  vel,  int   neval,
                       float *scale,  float *interv); 



void main(int argc, char *argv[])
/*===========================================================================*
 * Driver.                                                                   *
 *===========================================================================*/
{ 
  float  Fstart = 0.0,
         Fend   = 0.0,
         Fsamp  = 1.0,
         Diam   = 1.0,
         Vel    = 1.0;
  
  int    NPoles = DEFPOL;
  float *Poles;
  float  NormFactor, Gain, FreqInt, Freq, FDelta, SDen;

  int    NData;
  float *Data = fvector(1, MAXPNT);

  int    i, NEval = DEFEVL;
  
  FILE  *fp_in = stdin;

  
  getargs(argc, argv, &fp_in, 
	  &Fstart, &Fend, &Fsamp, &Diam, &Vel, &NEval, &NPoles);

  Poles = fvector(1, NPoles);

  doscale(Fstart, &Fend, Fsamp, Diam, Vel, NEval, &NormFactor, &FreqInt);

  getdata(fp_in, Data, &NData, MAXPNT);

  memcof(Data, NData, NPoles, &Gain, Poles);

  printhead(stdout, Fstart, Fend, NData, NPoles);
  
  for (i=0; i<NEval; i++) {
    Freq   = Fstart + i * FreqInt;
    FDelta = NormFactor * Freq;
    SDen   = 2.0 * NormFactor * evlmem(FDelta, Poles, NPoles, Gain);
    printf("%10.5g %10.5g\n", Freq, SDen);
  }
 
  exit(0);
}





static void doscale(float  fstart, float *fend, float fsamp,
		    float  dia,    float  vel,  int   neval,
                    float *scale,  float *interv )
/*===========================================================================*
 * Compute the scaling factor for frequencies and spectral densities, also   *
 * the frequency interval.                                                   *
 *===========================================================================*/
{
  *scale  = vel / (fsamp * dia);
  if (*fend <= 0.0)
    *fend = 0.5 / *scale;
  *interv = (*fend - fstart) / (neval - 1);

  if ((fstart < 0.0) || (*fend > fsamp / (2.0 * *scale)))
    message("memspec: doscale()", "frequency limit out of 0--Nyquist", ERROR);
}





static void getargs (int argc, char *argv[], FILE **fp_in, 
		     float *Fstart, float *Fend, float *Fsamp,
		     float *Diam,   float *Vel,  int *Neval, int *Npoles)
/*===========================================================================*
 * Do the command line.                                                      *
 *===========================================================================*/
{
  static char usage[] = 
    "usage: memspec [options] [input]\n\n"
    "where [options] are:\n"
    "-m #poles ... number of poles to fit                 [Default:  10]\n"
    "-n #evals ... number of points for evaluation        [Default: 512]\n"
    "-r fsamp  ... sampling frequency                     [Default: 1.0]\n"
    "-U vel    ... velocity scale for dimensionless freqs [Default: 1.0]\n"
    "-D dia    ... length   scale for dimensionless freqs [Default: 1.0]\n"
    "-s start  ... min frequency for spectral evaluation  [Default: 0.0]\n"
    "-e end    ... max frequency for spectral evaluation  [Def: Nyquist]\n";
  char   c;


  while (--argc && (*++argv)[0] == '-') {
    switch (c = *++argv[0]) {
    case 'h':
      fputs(usage, stderr);
      exit(0);
      break;
    case 'm':
      if (*++argv[0])
	*Npoles = atoi(*argv);
      else {
	*Npoles = atoi(*++argv);
	argc--;
      }
      break;
    case 'n':
      if (*++argv[0])
	*Neval = atoi(*argv);
      else {
	*Neval = atoi(*++argv);
	argc--;
      }
      break;
    case 'r':
      if (*++argv[0])
	*Fsamp = atof(*argv);
      else {
	*Fsamp = atof(*++argv);
	argc--;
      }
      break;
    case 'U':
      if (*++argv[0])
	*Vel = atof(*argv);
      else {
	*Vel = atof(*++argv);
	argc--;
      }
      break;
    case 'D':
      if (*++argv[0])
	*Diam = atof(*argv);
      else {
	*Diam = atof(*++argv);
	argc--;
      }
      break;
    case 's':
      if (*++argv[0])
	*Fstart = atof(*argv);
      else {
	*Fstart = atof(*++argv);
	argc--;
      }
      break;
    case 'e':
      if (*++argv[0])
	*Fend = atof(*argv);
      else {
	*Fend = atof(*++argv);
	argc--;
      }
      break;

    default:
      fprintf(stderr, "memspec: unknown option -- %c\n", c);
      exit(1);
      break;
    }
  }

  if (argc == 1)
    if ((*fp_in = fopen(*argv, "r")) == (FILE*) NULL) {
      fprintf(stderr, "memspec: couldn't open input file -- %s\n", *argv);
      exit(1);
    }
}




   
static void getdata (FILE *fp, float *Data, int *Ndata, int MAX)
/*===========================================================================*
 * Read data from file up to MAX points.                                     *
 *===========================================================================*/
{
  *Ndata = 0;
  while ((!feof(fp)) && *Ndata<MAX){
    ++*Ndata;
    fscanf(fp, "%f", Data + *Ndata);
  }
}





static void printhead(FILE *fp, float Fstart, float Fend, int Nd, int Np)
/*===========================================================================*
 * Print header statistics.                                                  *
 *===========================================================================*/
{

  static char *hdr_fmt[] = {
    "%-25s "      "Maximum entropy spectral analysis\n",
    "%-25s "      "Created\n",
    "%-25d "      "Number of data\n",
    "%-25d "      "Number of poles\n",
    "%-25.6g "    "Low frequency limit\n",
    "%-25.6g "    "High frequency limit\n"
    };
  static char name[] = "memspec";
  time_t tp;
  char   buf[BUFSIZ];


  tp = time((time_t*) NULL);
  strftime(buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));

  fprintf(fp, hdr_fmt[0], name);
  fprintf(fp, hdr_fmt[1], buf);
  fprintf(fp, hdr_fmt[2], Nd);
  fprintf(fp, hdr_fmt[3], Np);
  fprintf(fp, hdr_fmt[4], Fstart);
  fprintf(fp, hdr_fmt[5], Fend);
  fprintf(fp, "\n");
}




static float evlmem(float fdt, float d[], int m, float xms)
/*===========================================================================*
 * Given d[1..m], m, xms as returned by memcof(), this function returns the  *
 * power spectrum P(f) as a function of fdt = fDelta.                        *
 *===========================================================================*/
{
  int    i;
  float  sumr = 1.0, sumi = 0.0;
  double wr   = 1.0, wi   = 0.0,
         wpr, wpi, wtemp, theta;
    
  theta = TWOPI * fdt;
  wpr   = cos(theta);
  wpi   = sin(theta);

  for (i=1; i<=m; i++) {
    wr = (wtemp = wr)*wpr -    wi*wpi;
    wi =           wi*wpr + wtemp*wpi;
    sumr -= d[i]*wr;
    sumi -= d[i]*wi;
  }

  return xms/(sumr*sumr + sumi*sumi);
}



    
static void memcof(float data[], int n, int m, float *xms, float d[])
/*===========================================================================*
 * Given a real vector of data[1..n], and given m, this routine returns m    *
 * linear prediction coefficients as d[1..m], and returns the mean square    *
 * discrepancy as xms.                                                       *
 *===========================================================================*/
{
  int   k,j,i;
  float p = 0.0, *wk1, *wk2, *wkm;

    
  wk1 = fvector(1,n);
  wk2 = fvector(1,n);
  wkm = fvector(1,m);

  for (j=1; j<=n; j++) p += SQR(data[j]);
  *xms = p/n;

  wk1[1]   = data[1];
  wk2[n-1] = data[n];
  for (j=2; j<n; j++) {
    wk1[j]   = data[j];
    wk2[j-1] = data[j];
  }

  for (k=1; k<=m; k++) {
    float num = 0.0, denom = 0.0;

    for (j=1; j<=(n-k); j++) {
      num   += wk1[j] * wk2[j];
      denom += SQR(wk1[j]) + SQR(wk2[j]);
    }

    d[k]  = 2.0 * num/denom;
    *xms *= (1.0-SQR(d[k]));

    for (i=1; i<=(k-1); i++)
      d[i] = wkm[i] - d[k]*wkm[k-i];

    if (k == m) {
      free_fvector(wkm,1);
      free_fvector(wk2,1);
      free_fvector(wk1,1);
      return;
    }

    for (i=1; i<=k; i++)
      wkm[i] = d[i];
    for (j=1; j<=(n-k-1); j++) {
      wk1[j] -= wkm[k]*wk2[j];
      wk2[j]  = wk2[j+1] - wkm[k]*wk1[j+1];
    }
  }
}
