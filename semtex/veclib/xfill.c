/*****************************************************************************
 * xfill:   x[i] = alpha.                                                    *
 *****************************************************************************/


void dfill(int n, double alpha, double *x, int incx)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n ; i++ ) {
    *x = alpha;
    x += incx;
  }
}





void ifill(int n, int alpha, int *x, int incx)
{
  register int i;


  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n ; i++ ) {
    *x = alpha;
    x += incx;
  }
}





void sfill(int n, float alpha, float *x, int incx)
{
  register int i;


  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n ; i++ ) {
    *x = alpha;
    x += incx;
  }
}
