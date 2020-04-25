/*****************************************************************************
 * xgathr:  vector gather:  z[i] = x[y[i]].
 *
 * $Id: xgathr.c,v 9.1 2019/05/30 06:36:13 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dgathr (int_t n, const double* x, const int_t* y, double* z)
{
  register int_t i;
  
  for (i = 0; i < n; i++) z[i] = x[y[i]];
}


void igathr (int_t n, const int_t* x, const int_t* y, int_t* z)
{
  register int_t i;
  
  for (i = 0; i < n; i++) z[i] = x[y[i]];
}


void sgathr (int_t n, const float* x, const int_t* y, float* z)
{
  register int_t i;
  
  for (i = 0; i < n; i++) z[i] = x[y[i]];
}
