/*
 * index of minimum value
 */

int idmin (int n, double *x, int incx)
{
  int i;
  int indx = ( n > 0 ) ? 0 : -1;
  double xmin = *x;

  for(i = 0;i < n;i++) {
    if( *x < xmin ) {
      xmin = *x;
      indx = i;
    }
    x += incx;
  }
  return indx;
}

int ismin (int n, float *x, int incx)
{
  int i;
  int indx = ( n > 0 ) ? 0 : -1;
  float xmin = *x;

  for(i = 0; i < n; i++) {
    if( *x < xmin ) {
      xmin = *x;
      indx = i;
    }
    x += incx;
  }
  return indx;
}

