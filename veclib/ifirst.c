/*****************************************************************************
 * ifirst:  index of first non-zero value in x.                              *
 *****************************************************************************/


int ifirst(int n, const int *x, int incx)
{ 
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    if (*x) return i;
    x += incx;
  }
  
  return -1;
}
