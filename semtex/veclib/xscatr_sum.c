/*****************************************************************************
 * xscatr_sum:  vector scatter with summation:  z[y[i]] += x[i].
 *
 * $Id$
 *****************************************************************************/

  
void dscatr_sum (int n, const double* x, const int* y, double* z)
{
  register int  i;

  for (i = 0; i < n; i++) {
    *(z + *y) += *x;
    x++;
    y++;
  }
}


void iscatr_sum (int n, const int* x, const int* y, int* z)
{
  register int  i;

  for (i = 0; i < n; i++) {
    *(z + *y) += *x;
    x++;
    y++;
  }
}


void sscatr_sum (int n, const float* x, const int* y, float* z)
{
  register int  i;

  for (i = 0; i < n; i++) {
    *(z + *y) += *x;
    x++;
    y++;
  }
}
