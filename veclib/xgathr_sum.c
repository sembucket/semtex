/*****************************************************************************
 *  xgathr_sum:  vector gather, with summation:  z[i] += x[y[i]].
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dgathr_sum (int n, const double* x, const int* y, double* z)
{
  register int  i;

  for (i = 0; i < n; i++) z[i] += x[y[i]];
}


void igathr_sum (int n, const int* x, const int* y, int* z)
{
  register int  i;

  for (i = 0; i < n; i++) z[i] += x[y[i]];
}


void sgathr_sum (int n, const float* x, const int* y, float* z)
{
  register int  i;

  for (i = 0; i < n; i++) z[i] += x[y[i]];
}
