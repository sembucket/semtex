/*
 * mxm() - Matrix Multiply (double precision)
 *
 * The following function computes the matrix product C = A * B.  The matrices
 * are assumed to be in row-major order and must occupy consecutive memory lo-
 * cations.  The input quantities are as follows:
 *
 *     mxm(A,nra,B,nca,C,ncb) 
 *
 *     A    ... double* ... source vector one
 *     nra  ... int     ... number of rows in A
 *     B    ... double* ... source operand two
 *     nca  ... int     ... number of columns in A
 *     C    ... double* ... result vector
 *     ncb  ... int     ... number of columns in B
 *
 * A more general matrix multiply with arbitrary skips is given in mxma().
 */

#include "veclib/veclib.h"

#if !defined(mxm) 

void mxm (double* A, int nra, double* B, int nca, double *C, int ncb)
{
  double *a = A;
  double *b = B;
  double *c = C;
  double sum;

  int i, j, k;

  for (i = 0; i < nra; i++) {
    for (j = 0; j < ncb; j++) {

      b   = B + j;                  /* Next column of B    */
      sum = 0.;                     /* Clear sum           */

      for(k = 0; k < nca; k++) {    /* ------------------- */
	sum += a[k] * (*b);         /* Inner product loop  */
	b   += ncb;                 /* ------------------- */
      }
      *c++   = sum;                 /* Store and increment */
    }
    a += nca;                       /* Next row of A       */
  }

  return;                           /*      All Done       */
}

#endif
