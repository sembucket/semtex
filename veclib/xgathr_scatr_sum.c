/*****************************************************************************
 * xgathr_scatr_sum: vector gather-scatter with summation: z[y[i]] += w[x[i]].
 *
 * NB:  It is assumed that this operation is vectorizable, i.e. that there
 * are no repeated indices in the indirection vector y.
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dgathr_scatr_sum (int_t n,
		       const double*  w,
		       const int_t* x, const int_t* y,
		             double*  z)
{
  register int_t i;

#if defined(__uxp__)
#pragma loop novrec z
#elif defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) z[y[i]] += w[x[i]];
}


void igathr_scatr_sum (int_t n,
		       const int_t* w,
		       const int_t* x, const int_t* y,
		             int_t* z)
{
  register int_t i;

#if defined(__uxp__)
#pragma loop novrec z
#elif defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) z[y[i]] += w[x[i]];
}


void sgathr_scatr_sum (int_t n,
		       const float*   w,
		       const int_t* x, const int_t* y,
		             float*   z)
{
  register int_t i;

#if defined(__uxp__)
#pragma loop novrec z
#elif defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) z[y[i]] += w[x[i]];
}
