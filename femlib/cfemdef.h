#ifndef FEMDEF_H
#define FEMDEF_H
/*****************************************************************************
 * Common names and typedefs for finite-element codes.
 *
 * $Id$
 *****************************************************************************/

#if defined(_SX)		/* NEC SX-4.        */
  typedef long int integer;
  typedef double   real;
#else                           /* Everything else. */
  typedef int     integer;
  typedef double  real;
#endif

typedef integer int_t;
typedef real    real_t;

typedef struct { real x, y, z; } Point;
typedef Point                    Vector;

#define STR_MAX    256
#define F77NAME(x) x##_

typedef enum quadrature_kind {
  GJ  = 'G',	/* Gauss-Jacobi quadrature.         */
  GLJ = 'L',	/* Gauss-Lobatto-Jacobi quadrature. */
  GRJ = 'R',	/* Gauss-Radau-Jacobi quadrature.   */
  TRZ = 'U'	/* Trapezoidal (uniform mesh).      */
} QuadratureKind;

typedef enum solver_kind {
  DIRECT,	/* Cholesky back-substitution.                           */
  JACPCG,	/* Conjugate gradient, Jacobi (diagonal) preconditioner. */
  BLJPCG	/* Conjugate gradient, block Jacobi      preconditioner. */
} SolverKind;

typedef enum transform_kind {
  INVERSE = -1,	/* Basis    --> physical space. */
  FORWARD = +1	/* Physical --> basis    space. */
} TransformKind;

typedef enum exchange_kind {
  FULL,		/* Full exchange (preserves planar and element structure). */
  HALF		/* Half memory exchange across processes (planar only).    */
} ExchangeKind;

#endif