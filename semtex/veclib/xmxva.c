/*****************************************************************************
 * xmxva() - Matrix - Vector Multiply w/skips.                               *
 *                                                                           *
 * This following function computes the matrix-vector product C = A * B.     * 
 *                                                                           *
 *      mxva(A,iac,iar,B,ib,C,ic,nra,nca)                                    *
 *                                                                           *
 *      A   ... double* ... matrix factor (input)                            *
 *      iac ... int     ... increment in A between column elements           *
 *      iar ... int     ... increment in A between row elements              *
 *      B   ... double* ... vector factor (input)                            *
 *      ib  ... int     ... increment in B between consecutive elements      *
 *      C   ... double* ... vector product (output)                          *
 *      ic  ... int     ... increment in C between consecutive elements      *
 *      nra ... int     ... number of rows in A                              *
 *      nca ... int     ... number of columns in A                           *
 *                                                                           *
 * Consider BLAS2 xgemv as alternatives.                                     *
 *****************************************************************************/

#include <stdio.h>
#include <alplib.h>

#if !defined(mxva)


void dmxva(double *A, int iac, int iar, double *B,   int ib,
	   double *C, int ic,  int nra, int nca)
{
  register double *a, *b,
                  *c = C;
  register double  sum;
  register int     i, j;


  for(i = 0; i < nra; ++i) {
    sum = 0.0;
    a   =  A;
    A  += iac;
    b   =  B;
    for(j = 0; j < nca; ++j) {
      sum += (*a) * (*b); 
      a   += iar;
      b   += ib;
    }

    *c  = sum; 
     c += ic;
  }
}





void smxva(float *A, int iac, int iar, float *B,   int ib,
	   float *C, int ic,  int nra, int    nca)
{
  register float *a, *b,
                 *c = C;
  register float  sum;
  register int    i, j;

  for(i = 0; i < nra; ++i) {
    sum = 0.0F;
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
}

#endif
