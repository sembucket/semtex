/*****************************************************************************
 * xvpoly: Hoerner polynomial evaluation of a vector.
 *
 * m is the order of the polynomial and its coefficients are stored in c in
 * descending order: c[0] is the constant term.
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvpoly (integer n,
	     const double* x, integer incx, integer m,
	     const double* c, integer incc, 
	           double* y, integer incy)
{
  register integer i, j;
  register double  sum, xval;
  const    double  *csave, *cp;

  csave = c;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) {
    c    = (incc<0) ? csave+(-n+1)*incc : csave;
    cp   = c + incc;
    sum  = c[0];
    xval = x[i*incx];
    for (j = 0; j < m; j++) sum = sum * xval + cp[j*incc];
    y[i*incy] = sum;
  }
}


void svpoly (integer n,
	     const float* x, integer incx, integer m,
	     const float* c, integer incc, 
	           float* y, integer incy)
{
  register integer i, j;
  register float   sum, xval;
  const    float   *csave, *cp;

  csave = c;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) {
    c    = (incc<0) ? csave+(-n+1)*incc : csave;
    cp   = c + incc;
    sum  = c[0];
    xval = x[i*incx];
    for (j = 0; j < m; j++) sum = sum * xval + cp[j*incc];
    y[i*incy] = sum;
  }
}
