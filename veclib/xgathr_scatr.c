/*****************************************************************************
 * xgathr_scatr:  vector gather-scatter:  z[y[i]] = w[x[i]].
 *
 * NB: it is assumed that this operation is vectorizable, i.e. that there
 * are no repeated indices in the indirection vector y[i].
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif

  
void dgathr_scatr(int n, const double *w, const int *x,const int *y, double *z)
{
  register int i;

  for (i = 0; i < n; i++) z[y[i]] = w[x[i]];
}


void igathr_scatr(int n, const int *w, const int *x, const int *y, int *z)
{
  register int i;

  for (i = 0; i < n; i++) z[y[i]] = w[x[i]];
}


void sgathr_scatr(int n, const float *w, const int *x, const int *y, float *z)
{
  register int i;

  for (i = 0; i < n; i++) z[y[i]] = w[x[i]];
}
