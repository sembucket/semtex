/*****************************************************************************
 * xseq:  y[i] = alpha == x[i].                                              *
 *****************************************************************************/


void iseq(int n, int alpha, const int *x, int incx,
	                          int *y, int incy)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = (*x==alpha) ? 1 : 0;
    x += incx;
    y += incy;
  }
}
