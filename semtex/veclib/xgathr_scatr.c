/*****************************************************************************
 * xgathr_scatr:  vector gather-scatter:  z[y[i]] = w[x[i]].
 *
 * NB: it is assumed that this operation is vectorizable, i.e. that there
 * are no repeated indices in the indirection vector y[i].
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>

  
void dgathr_scatr (integer n,
		   const double* w,
		   const integer* x, const integer* y,
		   double* z)
{
  register integer i;

  for (i = 0; i < n; i++) z[y[i]] = w[x[i]];
}


void igathr_scatr (integer n,
		   const integer* w,
		   const integer* x, const integer* y,
		   integer* z)
{
  register integer i;

  for (i = 0; i < n; i++) z[y[i]] = w[x[i]];
}


void sgathr_scatr (integer n,
		   const float* w,
		   const integer* x, const integer* y,
		   float* z)
{
  register integer i;

  for (i = 0; i < n; i++) z[y[i]] = w[x[i]];
}
