/*****************************************************************************
 * dft_chk.c: check real-complex FFT, demonstrate scaling properties.
 * 
 * The Discrete Fourier Transform for complex data is taken as
 *
 *          1  N-1
 *   H(n) = -  Sum  h(k) exp (-i2PIkn/N),
 *          N  k=0
 *
 * with inverse
 *
 *             N-1
 *   h(k) =    Sum  H(n) exp (+i2PIkn/N).
 *             n=0
 *
 * The properties of the forward transform are such that the zero-frequency
 * value is twice the mean value of the data.  The lowest-frequency cosine
 * wave then has a value of +1.0 in the first real frequency location,
 * (i.e. data[1].Re = 1.0), while the lowest-frequency sine wave has a value
 * of -1.0 in the first imaginary frequency location (data[1].Im = -1.0).
 *
 * The normalizing factor that must be applied to the data on forward
 * transform is 2.0 / npts (or 1.0 / ncomplex).
 *
 * Note that for real data, the (real) Nyquist-frequency frequency point
 * is stored in the imaginary part of the zero-frequency transform location.
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "FFTutil.h"


static char prog[] = "dft_chk";


int main ()
/* ------------------------------------------------------------------------- *
 * Read number of data (power of two) on first line.
 * ------------------------------------------------------------------------- */
{ 
  int       i, m, npts;
  complex*  data;
  complex*  Wtab;
  real      norm;

  scanf ("%d", &npts);
  if (!ispow2 (npts)) message (prog, "input not a power of 2", ERROR);

  data = cvector (0, (npts/2)-1);
  Wtab = cvector (0, (npts/2)-1);

  for (i = 0; i < npts/2; i++) {
    scanf ((sizeof (real) == sizeof (float)) ? "%f" : "%lf", &data[i].Re);
    scanf ((sizeof (real) == sizeof (float)) ? "%f" : "%lf", &data[i].Im);
  }

  preFFT (Wtab, npts/2, -1);
  rcFFT  (data, npts/2, Wtab, npts/2, 1);

  printf ("-- Tranformed data:\n");
  norm = 2.0 / (npts);
  for (i = 0; i < npts/2; i++) {
    data[i].Re *= norm;
    data[i].Im *= norm;
    printf ("%f\t%f\n", data[i].Re, data[i].Im);
  }

  printf ("-- Inverted data:\n");
  rcFFT  (data, npts/2, Wtab, npts/2, 0);
  for (i = 0; i < npts/2; i++)
    printf ("%f\n%f\n", data[i].Re, data[i].Im);

  return EXIT_SUCCESS;
}
