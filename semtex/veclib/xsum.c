/*****************************************************************************
 * xsum:  sum = 0;  sum += x[i];                                             *
 *****************************************************************************/


double  dsum(int n, const double *x, int incx)
{
  register int     i;
  register double  sum;


  x += (incx<0) ? (-n+1)*incx : 0;
  sum = 0.0;

  for (i=0; i<n; i++) {
    sum += *x;
    x += incx;
  }
  
  return sum;
}





int  isum(int n, const int *x, int incx)
{
  register int  i, sum;


  x += (incx<0) ? (-n+1)*incx : 0;
  sum = 0;

  for (i=0; i<n; i++) {
    sum += *x;
    x += incx;
  }
  
  return sum;
}





float  ssum(int n, const float *x, int incx)
{
  register int    i;
  register float  sum;


  x += (incx<0) ? (-n+1)*incx : 0;
  sum = 0.0;

  for (i=0; i<n; i++) {
    sum += *x;
    x += incx;
  }
  
  return sum;
}
