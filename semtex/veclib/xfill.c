/*****************************************************************************
 * xfill:   x[i] = alpha.
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dfill (int n, double alpha, double *x, int incx)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha;
}


void ifill (int n, int alpha, int *x, int incx)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha;
}


void sfill (int n, float alpha, float *x, int incx)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha;
}
