/*****************************************************************************
 * xvpow:  y[i] = pow(x[i], alpha).
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dspow (const integer n, const double alpha,
	    const double* x, integer incx,
	          double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  
  for (i = 0; i < n; i++) y[i*incy] = pow (x[i*incx], alpha);
}


void sspow (const integer n, const float alpha,
	    const float* x, integer incx, 
	          float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  
  for (i = 0; i < n; i++)
#if defined(__GNUC__) || defined(__uxp__) || defined(_SX)
    y[i*incy] = (float) pow  (x[i*incx], alpha);
#else
    y[i*incy] =         powf (x[i*incx], alpha);
#endif
}
