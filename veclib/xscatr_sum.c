/*****************************************************************************
 * xscatr_sum:  vector scatter with summation:  z[y[i]] += x[i].
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

  
void dscatr_sum (int n, const double* x, const int* y, double* z)
{
  register int i;

  for (i = 0; i < n; i++) z[y[i]] += x[i];
}


void iscatr_sum (int n, const int* x, const int* y, int* z)
{
  register int i;

  for (i = 0; i < n; i++) z[y[i]] += x[i];
}


void sscatr_sum (int n, const float* x, const int* y, float* z)
{
  register int i;

  for (i = 0; i < n; i++) z[y[i]] += x[i];
}
