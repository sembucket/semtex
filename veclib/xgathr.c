/*****************************************************************************
 *  xgathr:  vector gather:  z[i] = x[y[i]].                                 *
 *****************************************************************************/


void dgathr (int n, const double *x, const int *y, double *z)
{
  register int  i;
  for (i = 0; i < n; i++) {*z = *(x + *y); y++; z++;}
}





void igathr (int n, const int *x, const int *y, int *z)
{
  register int  i;
  for (i = 0; i < n; i++) {*z = *(x + *y); y++; z++;}
}





void sgathr (int n, const float *x, const int *y, float *z)
{
  register int  i;
  for (i = 0; i < n; i++) {*z = *(x + *y); y++; z++;}
}
