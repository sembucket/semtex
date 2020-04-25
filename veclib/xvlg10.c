/*****************************************************************************
 * xvlg10:  y[i] = log10(x[i]).
 *
 * $Id: xvlg10.c,v 9.1 2019/05/30 06:36:14 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvlg10 (int_t n, const double* x, int_t incx,
                              double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = log10 (x[i*incx]);
}


void svlg10 (int_t n, const float* x, int_t incx,
                              float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

#if  defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) y[i*incy] = (float) log10  (x[i*incx]);
#else
  for (i = 0; i < n; i++) y[i*incy] =         log10f (x[i*incx]);
#endif
}
