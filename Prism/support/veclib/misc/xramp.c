/*
 * Ramp function
 */

void dramp (int n, double *ap, double *bp, double *x, const int incx)
{
  const double alpha = *ap;
  const double beta  = *bp;
  int i;
  for (i = 0; i < n; i++) {
    *x = alpha + (double) i * beta;
    x += incx;
  }
}

void sramp (int n, float *ap, float *bp, float *x, const int incx)
{
  const float alpha = *ap;
  const float beta  = *bp;
  int i;
  for (i = 0; i < n; i++) {
    *x = alpha + (float) i * beta;
    x += incx;
  }
}

void iramp (int n, int *ap, int *bp, int *x, const int incx)
{
  const double alpha = *ap;
  const double beta  = *bp;
  int i;
  for(i = 0; i < n; i++) {
    *x = alpha + i * beta;
    x += incx;
  }
}


