/*****************************************************************************
 * ixmax: index of maximum value in x.                                       *
 *****************************************************************************/


int idmax(int n, const double *x, int incx)
{
  register int     i, imax;
  register double  xmax;


  x += (incx<0) ? (-n+1)*incx : 0;
  xmax = *x;
  imax = 0;

  for (i=1; i<n; i++) {
    if (*x>xmax) {
      xmax = *x;
      imax = i;
    }
    x += incx;
  }

  return imax;
}





int iimax(int n, const int *x, int incx)
{
  register int  i, xmax, imax;


  x += (incx<0) ? (-n+1)*incx : 0;
  xmax = *x;
  imax = 0;

  for (i=1; i<n; i++) {
    if (*x>xmax) {
      xmax = *x;
      imax = i;
    }
    x += incx;
  }

  return imax;
}





int ismax(int n, const float *x, int incx)
{
  register int    i, imax;
  register float  xmax;


  x += (incx<0) ? (-n+1)*incx : 0;
  xmax = *x;
  imax = 0;

  for (i=1; i<n; i++) {
    if (*x>xmax) {
      xmax = *x;
      imax = i;
    }
    x += incx;
  }

  return imax;
}
