/*
 * Conjugate Gradient Solver
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"

/* Private functions */

static double *Jacobi (Field *U, BSystem *B);

/* ------------------------------------------------------------------------ *
 * Solve_CG() - Conjugate Gradient Solver                                   *
 *                                                                          *
 * This routine solves the algebraic system                                 *
 *                                                                          *
 *                               A u = F                                    *
 *                                                                          *
 * using a conjugate gradient iteration.  The input and output are the same *
 * as Solve().   The operation count is proportional to                     *
 *                                                                          *
 *                 CI  .=. J_eps [ K N^3 + K N^2 + K N ]                    *
 *                                                                          *
 * where         J_eps .=. sqrt (K N^3)                                     *
 *                                                                          *
 * is the number of iterations required for convergence.                    *
 * ------------------------------------------------------------------------ */

void Solve_CG (Field *U, Field *F, Bedge *Ubc, BSystem *B)
{
  double *M = Jacobi (U, B);

  Solve_CG_A (U, F, Ubc, B, Helmholtz, M);

  free (M);
  return;
}

/* General-purpose interface */

void Solve_CG_A 
    (Field *U, Field *F, Bedge *Ubc, BSystem *B, Operator A, double *M)
{
  double alpha, beta, eps, epso, epsr, tau;
  register int iter;

  const int nel  =  B->elements;
  const int bpts =  B->bpts;
  const int npts =  B->ipts + bpts;
  const int ntot =  U->nr * U->ns * nel;
  const int n    =  npts;

  /* Temporary arrays */

  tempVector (p, n);          /* Search direction     */
  tempVector (q, n);          /* Conditioned residual */
  tempVector (r, n);          /* Residual vector      */
  tempVector (v, n);          /* Solution             */

  tempVector (w, ntot);       /* Inner Product vector */

  /* ----------------------------------------------------------- *
   * Initialization                                              *
   *                                                             *
   * (1) Apply boundary conditions to the initial guess          *
   * (2) Compute the starting residual and mask it               *
   *                                                             *
   * ----------------------------------------------------------- */

  BC_set (Ubc, U, B);
  
  tau = Form_RHS (U, F, Ubc, B, A, r);
        Mask_RHS (U, B, r);

  dzero (n, v, 1);

  if (B->singular) dsadd (n, -dsum(n, r, 1)/n, r, 1, r, 1);

  iter = 0;
  epsr = sqrt(ddot (n, r, 1, r, 1));

  while (epsr > tau  && iter++ < n)
    {
      /* Preconditioner */

      dvmul (n, M, 1, r, 1, q, 1);


      /* Compute a new search direction */

      eps = ddot (n, r, 1, q, 1);
      if (iter == 1) {
	dcopy (n, q, 1, p, 1);
      } else {
	beta = eps / epso;
	dsvtvp (n, beta, p, 1, q, 1, p, 1);
      }
      

      /* ---------- Matrix Inner Product ---------- */

      Field_local (U, B, w, p);
      A           (U, B, w, w);
      Field_dsum  (U, B, w, q);
      Mask_RHS    (U, B, q);

      /* Update the solution and residual */

      alpha = eps / ddot(n, p, 1, q, 1);
      epso  = eps;

      daxpy (n, -alpha, q, 1, r, 1);
      daxpy (n,  alpha, p, 1, v, 1);

      epsr = sqrt(ddot(n, r, 1, r, 1));
    }

  /* Transform the correction to local format and add to U */

  Field_local (U, B, w, v);
  dvadd (ntot, w, 1, *U->field, 1, *U->field, 1);

  /* Free temporary vectors */

  freeVector (p);
  freeVector (q);
  freeVector (r);
  freeVector (v);
  freeVector (w); 
}


/*
 * Compute the Jacobi preconditioner for the Helmholtz equation
 */

static double *Jacobi (Field *U, BSystem *B)
{
  const int nr     = U->nr;
  const int ns     = U->ns;
  const int nrns   = U->nr * U->ns;
  const int nb     =(U->nr + U->ns - 2) << 1;
  const int nel    = B->elements;
  const int npts   = B->bpts + B->ipts;
  const int ni     = nrns - nb;
  double    lambda = B->constant;

  double **g1, **g2, **g3, **drt, **dst;
  int     *emap, *bmap;
  register int i, j, k;

  Matrix_IP *IP = B->CG;

  double *M = dvector (0, npts), *dM, *Mi;
  double  Mk [_MAX_NORDER * _MAX_NORDER],
          tmp[_MAX_NORDER];

  getops (nr, 0, 0, 0, &drt);
  getops (ns, 0, 0, 0, &dst);

  Mi = M + B->bpts;
  dzero (npts, M, 1);

  for (k = 0; k < nel; k++) {

    g1 = IP[k].g1;
    g2 = IP[k].g2;
    g3 = IP[k].g3;
    dM = Mk;

    dzero (nrns, dM, 1);

    for (i = 0; i < ns; i++) {
      for (j = 0; j < nr; j++, dM++) {
	dvmul(nr, drt[j], 1, drt[j], 1, tmp, 1);
	*dM    = ddot  (nr,  g1 [i], 1, tmp, 1);	  
	dvmul(ns, dst[i], 1, dst[i], 1, tmp, 1);
	*dM   += ddot  (ns,  *g2+j, nr, tmp, 1);

	if (g3)
	  *dM += 2. * g3[i][j] * drt[j][j] * dst[i][i];
	
	*dM += lambda * U[k].mass[i][j];
      }
    }
    
    emap = U[k].emap;
    bmap = B->bmap[k];
    for (i = 0; i < nb; i++)
      M [bmap[i]] += Mk [emap[i]];

    dgathr (ni, Mk, emap + nb, Mi);
    Mi += ni;
  }
  
  /* Invert and make sure it's positive-definite... */

  dvabs (npts, M, 1, M, 1);
  dvrecp(npts, M, 1, M, 1);

  return M;
}
