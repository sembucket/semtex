/*****************************************************************************
 * xgathr_scatr_sum: vector gather-scatter with summation: z[y[i]] += w[x[i]].
 *
 * NB:  It is assumed that this operation is vectorizable, i.e. that there
 * are no repeated indices in the indirection vector y.
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dgathr_scatr_sum (integer n,
		       const double*  w,
		       const integer* x, const integer* y,
		             double*  z)
{
  register integer i;

#if defined(__uxp__)
#pragma loop novrec z
#elif defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) z[y[i]] += w[x[i]];
}


void igathr_scatr_sum (integer n,
		       const integer* w,
		       const integer* x, const integer* y,
		             integer* z)
{
  register integer i;

#if defined(__uxp__)
#pragma loop novrec z
#elif defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) z[y[i]] += w[x[i]];
}


void sgathr_scatr_sum (integer n,
		       const float*   w,
		       const integer* x, const integer* y,
		             float*   z)
{
  register integer i;

#if defined(__uxp__)
#pragma loop novrec z
#elif defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) z[y[i]] += w[x[i]];
}
