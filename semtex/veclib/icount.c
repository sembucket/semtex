/*****************************************************************************
 * icount:  number of non-zero values in x.                                  *
 *****************************************************************************/

int icount(int n, const int *x, int incx)
{
  register int  i, sum=0;


  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++ ) {
    sum += (*x) ? 1 : 0;
    x += incx;
  }

  return sum;
}
