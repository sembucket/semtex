/*****************************************************************************
 * xvmul:  z[i] = x[i] * y[i].
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dvmul (int n, const double *x, int incx, const double *y, int incy,
	                                            double *z, int incz)
{
  register int  i;
  
  if (incx == 1 && incy == 1 && incz == 1)
    for (i = 0; i < n; i++) z[i] = x[i] * y[i];
  
  else {
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++) z[i*incz] = x[i*incx] * y[i*incy];
  }
}


void ivmul (int n, const int *x, int incx, const int *y, int incy,
	                                         int *z, int incz)
{
  register int  i;
  
  if (incx == 1 && incy == 1 && incz == 1)
    for (i = 0; i < n; i++) z[i] = x[i] * y[i];

  else {
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++) z[i*incz] = x[i*incx] * y[i*incy];
  }
}


void svmul (int n, const float *x, int incx, const float *y, int incy,
	                                           float *z, int incz)
{
  register int  i;

  if (incx == 1 && incy == 1 && incz == 1)
    for (i = 0; i < n; i++) z[i] = x[i] * y[i];

  else {
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++) z[i*incz] = x[i*incx] * y[i*incy];
  }
}
