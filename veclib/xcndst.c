/*****************************************************************************
 * xcndst:  conditional assignment:  if (y[i]) z[i] = x[i].                  *
 *****************************************************************************/


void dcndst(int n, const double *x, int incx, const int    *y, int incy,
                                                    double *z, int incz)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i=0; i<n; i++) {
    if (*y) *z = *x;
    x += incx;
    y += incy;
    z += incz;
  }
}





void icndst(int n, const int *x, int incx, const int  *y, int incy,
                                                 int  *z, int incz)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i=0; i<n; i++) {
    if (*y) *z = *x;
    x += incx;
    y += incy;
    z += incz;
  }
}





void scndst(int n, const float *x, int incx, const int   *y, int incy,
                                                   float *z, int incz)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i=0; i<n; i++) {
    if (*y) *z = *x;
    x += incx;
    y += incy;
    z += incz;
  }
}
