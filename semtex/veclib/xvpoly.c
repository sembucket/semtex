/*****************************************************************************
 * xvpoly: Hoerner polynomial evaluation of a vector.                        *
 *                                                                           *
 * m is the order of the polynomial and its coefficients are stored in c in  *
 * descending order: c[0] is the constant term.                              *
 *****************************************************************************/


void dvpoly(int n, const double *x, int incx, int m,
	           const double *c, int incc, 
		         double *y, int incy)
{
  register int    i,   j;
  register double sum, xval, *cp, *csave;


  csave = c;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    c    = (incc<0) ? csave+(-n+1)*incc : csave;
    cp   = c + incc;
    sum  = *c;
    xval = *x;
    for (j=0; j<m; j++) {
      sum = sum * xval + *cp;
      cp += incc;
    }
    *y = p;
    x += incx;
    y += incy;
  }
}





void svpoly(int n, const float  *x, int incx, int m,
	           const float  *c, int incc, 
		         float  *y, int incy)
{
  register int    i,   j;
  register float  sum, xval, *cp, *csave;


  csave = c;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    c    = (incc<0) ? csave+(-n+1)*incc : csave;
    cp   = c + incc;
    sum  = *c;
    xval = *x;
    for (j=0; j<m; j++) {
      sum = sum * xval + *cp;
      cp += incc;
    }
    *y = p;
    x += incx;
    y += incy;
  }
}
