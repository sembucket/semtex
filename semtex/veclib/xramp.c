/*****************************************************************************
 * xramp:  ramp function:  x[i] = alpha + i*beta.                            *
 *****************************************************************************/


void dramp(int n, double alpha, double beta, double *x, int incx)
{
  register int  i;

  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = alpha + i * beta;
    x += incx;
  }
}





void iramp(int n, int alpha, int beta, int *x, int incx)
{
  register int  i;

  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = alpha + i * beta;
    x += incx;
  }
}





void sramp(int n, float alpha, float beta, float *x, int incx)
{
  register int  i;

  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = alpha + i * beta;
    x += incx;
  }
}
