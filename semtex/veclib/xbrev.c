/*****************************************************************************
 *  xbrev:  byte-reversal routines.                                          *
 *****************************************************************************/


void dbrev (int n, const double *x, int incx, double *y, int incy)
{
  register char *cx, *cy, d;
  register int   i, j;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i<n; i++) {
    cx = (char*) x;
    cy = (char*) y;
    for (j = 0; j < 4; j++) { 
      d       = cx[j];
      cy[j]   = cx[7-j]; 
      cy[7-j] = d;
    }
    x += incx;
    y += incy;
  }
}





void ibrev (int n, const int *x, int incx, int *y, int incy)
{
  register char *cx, *cy, d;
  register int   i, j;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i<n; i++) {
    cx = (char*) x;
    cy = (char*) y;
    for (j = 0; j < 2; j++) { 
      d       = cx[j];
      cy[j]   = cx[3-j]; 
      cy[3-j] = d;
    }
    x += incx;
    y += incy;
  }
}





void sbrev (int n, const float *x, int incx, float *y, int incy)
{
  register char *cx, *cy, d;
  register int   i, j;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    cx = (char*) x;
    cy = (char*) y;
    for (j = 0; j < 2; j++) { 
      d       = cx[j];
      cy[j]   = cx[3-j]; 
      cy[3-j] = d;
    }
    x += incx;
    y += incy;
  }
}

