/*****************************************************************************
 * xscatr:  vector scatter:  z[y[i]] = x[i].                                 *
 *****************************************************************************/

  
void dscatr(int n, const double *x, const int *y, double *z)
{
  register int  i;


  for(i=0; i<n; i++) {
    *(z + *y) = *x;
    x++;
    y++;
  }
}





void iscatr(int n, const int *x, const int *y, int *z)
{
  register int  i;


  for(i=0; i<n; i++) {
    *(z + *y) = *x;
    x++;
    y++;
  }
}





void sscatr(int n, const float *x, const int *y, float *z)
{
  register int  i;


  for(i=0; i<n; i++) {
    *(z + *y) = *x;
    x++;
    y++;
  }
}
