/*****************************************************************************
 * xvpow:  y[i] = pow(x[i], alpha).
 *
 * $Id$
 *****************************************************************************/

#include <math.h>

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dspow (int n, const double alpha,
	    const double* x, int incx,
	    const double* y, int incy);
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  
  for (i = 0; i < n; i++) y[i*incy] = pow (x[i*incx], alpha);
}


void sspow (int n, const float alpha,
	    const float* x, int incx,
	    const float* y, int incy);
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  
  for (i = 0; i < n; i++)
#if (defined(__GNUC__) || defined (__uxp__)) /* -- No single-precision maths */
    y[i*incy] = (float) pow  (x[i*incx], alpha);
#else
    y[i*incy] =         powf (x[i*incx], alpha);
#endif
}
