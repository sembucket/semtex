/*
 * mxva() - Matrix - Vector Multiply w/skips (double precision)
 *
 * This following function computes the matrix-vector product C = A * B. 
 *
 *      mxva(A,iac,iar,B,ib,C,ic,nra,nca)
 *
 *      A   ... double* ... matrix factor (input)
 *      iac ... int     ... increment in A between column elements
 *      iar ... int     ... increment in A between row elements
 *      B   ... double* ... vector factor (input)
 *      ib  ... int     ... increment in B between consectuve elements
 *      C   ... double* ... vector product (output)
 *      ic  ... int     ... increment in C between consectuve elements
 *      nra ... int     ... number of rows in A
 *      nca ... int     ... number of columns in A
 *
 */

#include "veclib/veclib.h"

#if !defined(mxva)

void mxva(double* A,int iac,int iar,double* B,int ib,double *C,int ic,
	  int nra, int nca)
{
  register double *a, *b,
                  *c = C;
  register double  sum;
  register int     i, j;

  for(i = 0; i < nra; ++i) {
    sum = 0.;
    a   = A;
    A  += iac;
    b   = B;
    for(j = 0; j < nca; ++j) {
      sum += (*a) * (*b); 
      a   += iar;
      b   += ib;
    }

    *c  = sum; 
     c += ic;
  }

  return;
}

#endif
