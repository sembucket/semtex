/*****************************************************************************
 * lany: return 1 if any x are true: iany = 0; if (x[i]) lany = 1.           *
 *****************************************************************************/


int lany (int n, const int *x, int incx)
{ 
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) {
    if (*x) return 1;
    x += incx;
  }
  
  return 0;
}
