/*****************************************************************************
 * xvsgn:  y[i] = sign(x[i])  --  sign = -1 if x<0, else +1.
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvsgn (integer n, const double* x, integer incx,
                             double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (x[i*incx] < 0.0) ? -1.0 : 1.0;
}


void ivsgn (integer n, const integer* x, integer incx,
                             integer* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (x[i*incx] < 0) ? -1 : 1;
}


void svsgn (integer n, const float* x, integer incx,
                             float* y, integer incy)
{
  register integer i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (x[i*incx] < 0.0F) ? -1.0F : 1.0F;
}
