/*****************************************************************************
 * xscatr:  vector scatter:  z[y[i]] = x[i].                                 *
 *****************************************************************************/

  
void dscatr(int n, double *x, int *y, double *z)
{
  register int  i;


  for(i=0; i<n; i++) {
    *(z + *y) = *x;
    x++;
    y++;
  }
}





void iscatr(int n, int *x, int *y, int *z)
{
  register int  i;


  for(i=0; i<n; i++) {
    *(z + *y) = *x;
    x++;
    y++;
  }
}





void sscatr(int n, float *x, int *y, float *z)
{
  register int  i;


  for(i=0; i<n; i++) {
    *(z + *y) = *x;
    x++;
    y++;
  }
}
