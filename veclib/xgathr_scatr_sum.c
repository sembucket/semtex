/*****************************************************************************
 * xgathr_scatr_sum: vector gather-scatter with summation: z[y[i]] += w[x[i]].
 *
 * NB:  It is assumed that this operation is vectorizable, i.e. that there
 * are no repeated indices in the indirection vector y.
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif

  
void dgathr_scatr_sum (int n, const double* w, const int*    x, 
		              const int*    y,       double* z)
{
  register int i;

  for (i = 0; i < n; i++) z[y[i]] += w[x[i]];
}


void igathr_scatr_sum (int n, const int* w, const int* x,
		              const int* y,       int* z)
{
  register int i;

  for (i = 0; i < n; i++) z[y[i]] += w[x[i]];
}


void sgathr_scatr_sum (int n, const float* w, const int*   x,
		              const int*   y,       float* z)
{
  register int i;

  for (i = 0; i < n; i++) z[y[i]] += w[x[i]];
}
