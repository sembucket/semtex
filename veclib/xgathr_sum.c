/*****************************************************************************
 *  xgathr_sum:  vector gather, with summation:  z[i] += x[y[i]].
 *
 * $Id$
 *****************************************************************************/


void dgathr_sum (int n, const double* x, const int* y, double* z)
{
  register int  i;

  for (i = 0; i < n; i++) {
    *z += *(x + *y);
    y++;
    z++;
  }
}


void igathr_sum (int n, const int* x, const int* y, int* z)
{
  register int  i;

  for (i = 0; i < n; i++) {
    *z += *(x + *y);
    y++;
    z++;
  }
}


void sgathr_sum (int n, const float* x, const int* y, float* z)
{
  register int  i;

  for (i = 0; i < n; i++) {
    *z += *(x + *y);
    y++;
    z++;
  }
}
