/*****************************************************************************
 *  xgathr:  vector gather:  z[i] = x[y[i]].                                 *
 *****************************************************************************/


void dgathr(int n, double *x, int *y, double *z)
{
  register int  i;

  
  for (i=0; i<n; i++) {
    *z = *(x + *y);
    y++;
    z++;
  }
}





void igathr(int n, int *x, int *y, int *z)
{
  register int  i;

  
  for (i=0; i<n; i++) {
    *z = *(x + *y);
    y++;
    z++;
  }
}





void sgathr(int n, float *x, int *y, float *z)
{
  register int  i;

  
  for (i=0; i<n; i++) {
    *z = *(x + *y);
    y++;
    z++;
  }
}
