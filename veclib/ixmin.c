/*****************************************************************************
 * ixmin:  index of minimum value in x.                                      *
 *****************************************************************************/


int idmin(int n, const double *x, int incx)
{
  register int     i, imin;
  register double  xmin;


  x += (incx<0) ? (-n+1)*incx : 0;
  xmin = *x;
  imin = 0;

  for (i=1; i<n; i++) {
    if (*x<xmin) {
      xmin = *x;
      imin = i;
    }
    x += incx;
  }

  return imin;
}





int iimin(int n, const int *x, int incx)
{
  register int  i, xmin, imin;


  x += (incx<0) ? (-n+1)*incx : 0;
  xmin = *x;
  imin = 0;

  for (i=1; i<n; i++) {
    if (*x<xmin) {
      xmin = *x;
      imin = i;
    }
    x += incx;
  }

  return imin;
}





int ismin(int n, const float *x, int incx)
{
  register int    i, imin;
  register float  xmin;


  x += (incx<0) ? (-n+1)*incx : 0;
  xmin = *x;
  imin = 0;

  for (i=1; i<n; i++) {
    if (*x<xmin) {
      xmin = *x;
      imin = i;
    }
    x += incx;
  }

  return imin;
}
