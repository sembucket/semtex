/*****************************************************************************
 * xneg:  x[i] = -x[i].
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dneg (int n, double *x, int incx)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = -x[i*incx];
}


void ineg (int n, int *x, int incx)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = -x[i*incx];
}


void sneg (int n, float *x, int incx)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = -x[i*incx];
}
