/*****************************************************************************
 * xneg:  x[i] = -x[i].                                                      *
 *****************************************************************************/


void dneg(int n, double *x, int incx)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = -*x;
    x += incx;
  }
}





void ineg(int n, int *x, int incx)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = -*x;
    x += incx;
  }
}





void sneg(int n, float *x, int incx)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = -*x;
    x += incx;
  }
}
