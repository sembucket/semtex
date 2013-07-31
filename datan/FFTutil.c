/*****************************************************************************
 *           UTILITIES AND PROCEDURES FOR 1-DIMENSIONAL FFTS.
 *
 * See `Numerical Recipes in C', 2nd edn, by W.H. Press et al., CUP.
 *
 * $Id$
 *****************************************************************************/
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "FFTutil.h"


static complex W;
#define SWAP(a, b)  W = (a); (a) = (b); (b) = W


void message (const char *routine, const char *text, int level)
/* ------------------------------------------------------------------------- *
 * A simple error handler.
 * ------------------------------------------------------------------------- */
{
  switch (level) {
  case WARNING:
    fprintf (stderr, "WARNING: %s: %s\n", routine, text); 
    break;
  case ERROR:
    fprintf (stderr, "ERROR: %s: %s\n", routine, text); 
    break;
  case REMARK:
    fprintf (stdout, "%s: %s\n", routine, text);
    break;
  default:
    fprintf (stderr, "bad error level in message: %d\n", level);
    exit (EXIT_FAILURE);
    break;
  }

  if (level == ERROR) exit (EXIT_FAILURE);
}


complex* cvector (long nl, long nh)
/* ------------------------------------------------------------------------- *
 * Allocates a complex vector with subscript range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  complex *v;

  v = (complex*) malloc ((size_t) ((nh-nl+1) * sizeof (complex)));
  if (v) return v-nl;

  message ("cvector()", "allocation failure", WARNING);
  return NULL;
}


int* ivector (long nl, long nh)
/* ------------------------------------------------------------------------- *
 * Allocates an int vector with subscript range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  int *v;

  v = (int*) malloc ((size_t) ((nh-nl+1) * sizeof (int)));
  if (v) return v-nl;

  message ("ivector()", "allocation failure", WARNING);
  return NULL;
}


real* rvector (long nl, long nh)
/* ------------------------------------------------------------------------- *
 * Allocates a real vector with subscript range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  real  *v;

  v = (real*) malloc ((size_t) ((nh-nl+1) * sizeof (real)));
  if (v) return v-nl;

  message ("rvector()", "allocation failure", WARNING);
  return NULL;
}


void freeCvector (complex *v, long nl)
/* ------------------------------------------------------------------------- *
 * Frees a complex vector allocated by cvector().
 * ------------------------------------------------------------------------- */
{
  free (v + nl);
}


void freeRvector (real *v, long nl)
/* ------------------------------------------------------------------------- *
 * Frees a real vector allocated by rvector().
 * ------------------------------------------------------------------------- */
{
  free (v + nl);
}


void freeIvector (int *v, long nl)
/* ------------------------------------------------------------------------- *
 * Frees an int vector allocated by ivector().
 * ------------------------------------------------------------------------- */
{
  free (v + nl);
}


int ispow2 (int k)
/* ------------------------------------------------------------------------- *
 * Is k a strictly positive integer power of two?
 * ------------------------------------------------------------------------- */
{
  while ((k%2) == 0 && k>2) k /= 2;
  if (k != 2) return 0;
  else        return 1;
}


int roundpow2 (int k)
/* ------------------------------------------------------------------------- *
 * Return the power of two greater than or equal k.
 * ------------------------------------------------------------------------- */
{
  if (k > 1) while (!ispow2 (k)) k++;

  return k;
}

 
void preFFT (complex *Wtab, const int K, const int Sign)
/* ------------------------------------------------------------------------- *
 * Make angular factors for FFT in Wtab[0..N-1].
 * K is the length of FFT; i.e., the number of complex values.
 * Sign is -1 for exp(-ikx) on forward transform (the opposite to Numerical
 * Recipes; use 1 to get their definition).
 * ------------------------------------------------------------------------- */
{
  int     i;
  double  theta;

  Wtab[0].Re = 1.0;   Wtab[0].Im = 0.0;

  for (i = 1; i < K; i++) {
    theta = i * M_PI / (double) K;
    Wtab[i].Re =        cos(theta);
    Wtab[i].Im = Sign * sin(theta);
  }
}


void cFFT (complex*        Data   ,
	   const int       N      ,
	   const complex*  Wtab   ,
	   const int       TabLen ,
	   const int       Forward)
/* ------------------------------------------------------------------------- *
 * Complex-complex FFT.   Forward = 1 ==> forward FFT, 0 ==> inverse.
 *
 * The array Data is assumed to be indexed starting at zero.
 * Length N is the number of complex data in transform (power of 2).
 * No normalization done here; divide by N as appropriate.
 * Angular factors must be precomputed (by preFFT()).
 * ------------------------------------------------------------------------- */
{
  register int	 mmax, m, p, q, s, t, tstep;
  register real  tempr, tempi;
  const    int   Non2 = N >> 1, Nm = N - 1;
  
  s = 0;			/* -- Bit-reversal. */
  for (t =0; t < Nm; t++) {
    if (s > t) {
      SWAP(Data[s], Data[t]);
    }
    m = Non2;
    while (m <= s) {
      s -= m;
      m >>= 1;
    }
    s += m;
  }

  mmax = 1;
  q    = TabLen;
  while (N > mmax) {		/* -- Butterfly. */
    p = 0;
    tstep = mmax << 1;
    for (m = 0; m < mmax; m++) {
      for (t = m; t < N; t += tstep) {
	s = t + mmax;
	tempr = Wtab[p].Re * Data[s].Re;
	tempi = Wtab[p].Re * Data[s].Im;
	if (Forward) {
	  tempr -= Wtab[p].Im * Data[s].Im;
	  tempi += Wtab[p].Im * Data[s].Re;
	} else {
	  tempr += Wtab[p].Im * Data[s].Im;
	  tempi -= Wtab[p].Im * Data[s].Re;
	}  
	Data[s].Re = Data[t].Re - tempr;
	Data[s].Im = Data[t].Im - tempi;
	Data[t].Re += tempr;
	Data[t].Im += tempi;
      }
      p += q;
    }
    mmax = tstep;
    q >>= 1;
  }
}


void rcFFT (complex*        Data   ,
	    const int       N      ,
	    const complex*  Wtab   ,
	    const int       TabLen ,
	    const int       Forward)
/* ------------------------------------------------------------------------- *
 * Perform FFT of real-stored-as-complex data (Real-complex FFT).
 * Forward = 1 ==> forward FFT, 0 ==> inverse.
 *
 * The array Data is assumed to be indexed starting at zero, with N complex
 * data (power of 2), i.e., 2N real values in total.
 *
 * No normalization done; divide by N as appropriate.
 * ------------------------------------------------------------------------- */
{
  register int	 k, revk;
  register real  c2, h1r, h1i, h2r, h2i;
  const    int   Non2 = N >> 1;
  const    real  c1   = 0.5;
  
  if (Forward) {
    cFFT (Data, N, Wtab, TabLen, Forward);
    c2 = -c1;
  }
  else
    c2 = c1;

  for (k = 1; k < Non2; k++) {
    revk = N - k;
    
    h1r =  c1 * (Data[k].Re + Data[revk].Re);
    h1i =  c1 * (Data[k].Im - Data[revk].Im);
    h2r = -c2 * (Data[k].Im + Data[revk].Im);
    h2i =  c2 * (Data[k].Re - Data[revk].Re);
    
    Data[k].Re    =  h1r + Wtab[k].Re * h2r;
    Data[k].Im    =  h1i + Wtab[k].Re * h2i;
    Data[revk].Re =  h1r - Wtab[k].Re * h2r;
    Data[revk].Im = -h1i + Wtab[k].Re * h2i;
    if (Forward) {
      Data[k].Re    -= Wtab[k].Im * h2i;
      Data[k].Im    += Wtab[k].Im * h2r;
      Data[revk].Re += Wtab[k].Im * h2i;
      Data[revk].Im += Wtab[k].Im * h2r;
    } else {
      Data[k].Re    += Wtab[k].Im * h2i;
      Data[k].Im    -= Wtab[k].Im * h2r;
      Data[revk].Re -= Wtab[k].Im * h2i;
      Data[revk].Im -= Wtab[k].Im * h2r;
    }
  }
  Data[Non2].Im = -Data[Non2].Im;
  
  if (Forward) {
    Data[0].Re = (h1r=Data[0].Re) + Data[0].Im;
    Data[0].Im =  h1r             - Data[0].Im;
  }
  else {
    Data[0].Re = c1 * ((h1r = Data[0].Re) + Data[0].Im);
    Data[0].Im = c1 * ( h1r - Data[0].Im);
    cFFT (Data, N, Wtab, TabLen, 0);
  }
}

 
void pcFFT (complex*   Zbuf   ,
	    const int  N      ,
	    const int  Forward)
/* ------------------------------------------------------------------------- *
 * Real-complex FFT of two lots of real data.
 *
 * Parameter N is the number of complex data which would exist in each of
 * the two DFTs if they were not packed.  The +ve-frequency part of the
 * first DFT is packed into the low half of Zbuf, with the real values from
 * the Nyquist frequency stored in the imaginary part of the  zero-frequency
 * allocation.  The +ve-frequency part of the second DFT is stored in the
 * upper half of Zbuf, but reflected, so that the zero and Nyquist frequency
 * data are located at the high end of Zbuf, with progressively higher
 * frequencies towards the centre of Zbuf.
 *
 * To carry out the forward FFT of two lots of real data, the first lot of
 * data are stored in the real locations in Zbuf, while the second lot are
 * are stored in the imaginary locations.  Forward FFT Zbuf then call
 * pcFFT(Zbuf, N, 1).  Unscrambling of the transform is carried out until
 * the +ve-frequency parts of the DFTs are separated as described above.
 *
 * When using the inverse process, the two DFTs are assumed to be packed as
 * described.  First call pcFFT(Zbuf, N, 0), which does the mixing-up of the
 * two DFTs, then do inverse FFT.  The two lots of real data are in the real
 * and imaginary locations of Zbuf.
 *
 * References: Numerical Recipes sect 12.3, Bendat & Piersol 1971 sect 9.84.
 * ------------------------------------------------------------------------- */
{
  register int      k, revk;
  register real     s1, s2, s3;
  register complex  A, B;
  const    int      Non2 = N >> 1;

  if (Forward) {	/* -- Take mixed-up DFTs & unscramble. */
    s1 = Zbuf[Non2].Re;
    s2 = Zbuf[Non2].Im;
    s3 = Zbuf[   0].Im;

    for (k = Non2 - 1; k > 0; k--) {
      revk = N - k;

      A = Zbuf[k];
      B = Zbuf[revk];
      Zbuf[k].Re      = 0.5 * (A.Re + B.Re);
      Zbuf[k].Im      = 0.5 * (A.Im - B.Im);
      Zbuf[revk-1].Re = 0.5 * (A.Im + B.Im);
      Zbuf[revk-1].Im = 0.5 * (B.Re - A.Re);
    }

    Zbuf[0].Im   = s1;
    Zbuf[N-1].Re = s3;
    Zbuf[N-1].Im = s2;

  } else {			/* -- Take DFTs of real data & scramble up. */
    s1 = Zbuf[N-1].Re;
    s2 = Zbuf[N-1].Im;
    s3 = Zbuf[0].Im;

    for (k = 1; k < Non2; k++) {
      revk = N - k;

      A = Zbuf[k];
      B = Zbuf[revk-1];
      Zbuf[k].Re    = A.Re - B.Im;
      Zbuf[k].Im    = A.Im + B.Re;
      Zbuf[revk].Re = A.Re + B.Im;
      Zbuf[revk].Im = B.Re - A.Im;
    }
    
    Zbuf[0].Im    = s1;
    Zbuf[Non2].Re = s3;
    Zbuf[Non2].Im = s2;
  }
}
