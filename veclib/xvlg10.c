/*****************************************************************************
 * xvlg10:  y[i] = log10(x[i]).
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif


void dvlg10 (integer n, const double* x, integer incx,
                              double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = log10 (x[i*incx]);
}


void svlg10 (integer n, const float* x, integer incx,
                              float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

#if defined(__GNUC__) || defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) y[i*incy] = (float) log10  (x[i*incx]);
#else
  for (i = 0; i < n; i++) y[i*incy] =         log10f (x[i*incx]);
#endif
}
