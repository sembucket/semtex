/*****************************************************************************
 * xmxva() - Matrix - Vector Multiply w/skips.
 *
 * This following function computes the matrix-vector product C = A * B.
 *
 *      mxva(A,iac,iar,B,ib,C,ic,nra,nca)
 *
 *      A   ... double* ... matrix factor (input)
 *      iac ... int     ... increment in A between column elements
 *      iar ... int     ... increment in A between row elements
 *      B   ... double* ... vector factor (input)
 *      ib  ... int     ... increment in B between consecutive elements
 *      C   ... double* ... vector product (output)
 *      ic  ... int     ... increment in C between consecutive elements
 *      nra ... int     ... number of rows in A
 *      nca ... int     ... number of columns in A
 *
 * Consider BLAS2 xgemv as alternatives.
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <femdef.h>
#include <alplib.h>


void dmxva(double* A, integer iac, integer iar, double* B, integer ib,
	   double* C, integer ic,  integer nra, integer nca)
{
  register double  *a, *b,
                   *c = C;
  register double  sum;
  register integer i, j;


  for (i = 0; i < nra; ++i) {
    sum = 0.0;
    a   =  A;
    A  += iac;
    b   =  B;
    for (j = 0; j < nca; ++j) {
      sum += (*a) * (*b); 
      a   += iar;
      b   += ib;
    }

    *c  = sum; 
     c += ic;
  }
}


void smxva (float* A, integer iac, integer iar, float* B, integer ib,
	    float* C, integer ic,  integer nra, integer nca)
{
  register float   *a, *b,
                   *c = C;
  register float   sum;
  register integer i, j;

  for (i = 0; i < nra; ++i) {
    sum = 0.0F;
    a   = A;
    A  += iac;
    b   = B;
    for (j = 0; j < nca; ++j) {
      sum += (*a) * (*b); 
      a   += iar;
      b   += ib;
    }

    *c  = sum; 
     c += ic;
  }
}
