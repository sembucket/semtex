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

typedef enum basis_kind {
  STD = 'U',	/* Standard, uniformly-spaced mesh. */
  GLL = 'S'	/* Gauss-Lobatto-Legendre mesh.     */
} BasisKind;

typedef enum quadrature_kind {
  GL = 'G',	/* Gauss-Legendre quadrature.      */
  LL = 'L',	/* Lobatto-Legendre quadrature.    */
  GR = 'R'	/* Special Gauss-Radau quadrature. */
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
