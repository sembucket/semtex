/*****************************************************************************
 * xmxv: Matrix - Vector Multiply.                                           *
 *                                                                           *
 * This following function computes the matrix-vector product C = A * B.     *
 * The matrix A is assumed to be stored in row-major order and must occupy   *
 * consecutive memory locations.                                             *
 *                                                                           *
 *      mxv(A,nra,B,nca,C)                                                   *
 *                                                                           *
 *      A   ... double* ... matrix factor (input)                            *
 *      nra ... int     ... number of rows in A                              *
 *      B   ... double* ... vector factor (input)                            *
 *      nca ... int     ... number of columns in A                           *
 *      C   ... double* ... vector product (output)                          *
 *                                                                           *
 * A more general matrix-vector multiply with arbitrary skips is given in    *
 * mxva().  Consider also the BLAS2 routines xgemv.                          *
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dmxv (double* A, integer nra, double* B, integer nca, double* C)
{
  register double  *a = A,
                   *c = C;
  register double  sum;
  register integer i, j;

  for (i = 0; i < nra; i++) {
    sum  = 0.0;
    for (j = 0; j < nca; j++) sum += (*a++) * B[j]; 
    *c++ = sum;
  }
}


void smxv (float* A, integer nra, float* B, integer nca, float* C)
{
  register float   *a = A,
                   *c = C;
  register float   sum;
  register integer i, j;

  for (i = 0; i < nra; i++) {
    sum  = 0.0F;
    for (j = 0; j < nca; j++) sum += (*a++) * B[j]; 
    *c++ = sum;
  }
}
