/*****************************************************************************
 * xgathr_scatr:  vector gather-scatter:  z[y[i]] = w[x[i]].                 *
 *****************************************************************************/

  
void dgathr_scatr(int n, const double *w, const int *x,const int *y, double *z)
{
  register int  i;
  for (i = 0; i < n; i++) {*(z + *y) = *(w + *x); x++; y++;}
}





void igathr_scatr(int n, const int *w, const int *x, const int *y, int *z)
{
  register int i;
  for (i = 0; i < n; i++) {*(z + *y) = *(w + *x); x++; y++;}
}





void sgathr_scatr(int n, const float *w, const int *x, const int *y, float *z)
{
  register int  i;
  for (i = 0; i < n; i++) { *(z + *y) = *(w + *x); x++; y++; }
}
