/*
 * SOLVE_CG_Q - Conjugate Gradient Solver for the Helmholtz Equation
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
#include "mortar/mortar.h"

/* Private functions */

static double *Jacobi (Field *U, BSystem *M);
static void    dsavg  (Field *U, BSystem *M, double *field);


/* ------------------------------------------------------------------------ *
 * Solve_CG_q() - Conjugate Gradient Solver ... Mortar Element Version      *
 *                                                                          *
 * This routine solves the algebraic system corresponding to                *
 *                                                                          *
 *                   a(w,v) = -(w,f) - (w,h) + a(w,g)                       *
 *                                                                          *
 * using the conjugate gradient method.  The inputs are the initial guess   *
 * "g", the force "f", the boundary conditions and the global matrix        *
 * system.                                                                  *
 * ------------------------------------------------------------------------ */

void Solve_CG_q (Field *U, Field *F, Bedge *Ubc, BSystem *B) 
{
  double *M = Jacobi (U, B);

  Solve_CG_A_q  (U, F, Ubc, B, Helmholtz, M);

  free (M);
  return;
}

/* General-purpose interface */

void Solve_CG_A_q
    (Field *U, Field *F, Bedge *Ubc, BSystem *B, Operator A, double *M)
{
  double alpha, beta, eps, epso, epsr, tau;
  int iter;

  const int nel  =  B->elements;
  const int ntot =  U->nr * U->ns * nel;
  const int n    =  ntot;
  Patch    *P    =  (Patch*) B->other;

  /* Temporary arrays */

  tempVector (p, n);          /* Search direction     */
  tempVector (q, n);          /* Conditioned residual */
  tempVector (r, n);          /* Residual vector      */
  tempVector (v, n);          /* Solution             */

  tempVector (w, n);          /* Inner Product Vector */

  /* ----------------------------------------------------------- *
   * Initialization                                              *
   *                                                             *
   * (1) Apply boundary conditions to the initial guess          *
   * (2) Compute the starting residual and mask it               *
   *                                                             *
   * ----------------------------------------------------------- */

  BC_set (Ubc, U, B);
  
  tau = Form_RHS_q (U, F, Ubc, B, A, r);
        Mask_RHS_q (U, B, r);

  dzero (n, v, 1);

  if (B->singular) dsadd (n, -dsum(n,r,1)/n, r, 1, r, 1);

  iter = 0;
  epsr = sqrt(ddot(n, r, 1, r, 1));

  while (epsr > tau && iter++ < n)
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

      dcopy       (n, p, 1, w, 1);
      Project_m_s (U, P, w);                   /* w =         Z p   */
      A           (U, B, w, w);                /* w =       A Z p   */
      Project_s_m (U, P, w);                   /* w =    Z' A Z p   */
      dcopy       (n, w, 1, q, 1);             /*                   */
      Mask_RHS_q  (U, B, q);                   /* .... mask(q) .... */
      
      /* Update the solution and residual */

      alpha = eps / ddot(n, p, 1, q, 1);
      epso  = eps;

      daxpy (n, -alpha, q, 1, r, 1);
      daxpy (n,  alpha, p, 1, v, 1);

      epsr = sqrt(ddot(n, r, 1, r, 1));
    }
  
  Project_m_s (U, P, v);
  dvadd (n, v, 1, *U->field, 1, *U->field, 1);
  
  /* Free temporary vectors */

  freeVector (p);
  freeVector (q);
  freeVector (r);
  freeVector (v);
  freeVector (w);
}

/* -----------------  P R I V A T E    F U N C T I O N S  ----------------- */

/*
 * Compute the diagonal Helmholtz preconditioner
 */

static double *Jacobi (Field *U, BSystem *B)
{
  const int    nr     = U->nr,
           ns     = U->ns,
           nrns   = U->nr * U->ns,
           nel    = B->elements,
           ntot   = nrns * nel;
  double  *M      = dvector (0, ntot),
           lambda = B->constant;
  Patch   *P      = B->other;
  double **g1, **g2, **g3, **drt, **dst, *dp;
  register int i, j, k;

  Matrix_IP *IP = B->CG;

  double tmp [_MAX_NORDER];

  getops (nr, 0, 0, 0, &drt);
  getops (ns, 0, 0, 0, &dst);

  dp = M;

  for (k = 0; k < nel; k++) {

    g1 = IP[k].g1;
    g2 = IP[k].g2;
    g3 = IP[k].g3;
    
    for (i = 0; i < ns; i++) {
      for (j = 0; j < nr; j++, dp++) {
	dvmul(nr, drt[j], 1, drt[j], 1, tmp, 1);
	*dp    = ddot  (nr,  g1 [i], 1, tmp, 1);	  
	dvmul(ns, dst[i], 1, dst[i], 1, tmp, 1);
	*dp   += ddot  (ns,  *g2+j, nr, tmp, 1);

	if (g3)
	  *dp += 2. * g3[i][j] * drt[j][j] * dst[i][i];
	
	*dp += lambda * U[k].mass[i][j];
      }
    }
  }
  
  Project_s_m (U, P, M);     /* Project diagonal -> mortars  */
  Project_m_s (U, P, M);     /* Project back                 */
  dsavg       (U, B, M);     /* Average across boundaries    */

  /* Invert and make sure it's positive-definite... */

  dvabs (ntot, M, 1, M, 1);
  dvrecp(ntot, M, 1, M, 1);

  return M;
}

/*
 * Direct Stiffness Summation - Simple averaging
 */

static void dsavg (Field *U, BSystem *M, double *field)
{
  const int   bpts  = M->bpts,
          nel   = M->elements,
          ntot  = U->nr * U->ns,
          nb    =(U->nr + U->ns - 2) << 1;
  int**   bmap  = M->bmap;
  double  *fptr;
  register int i, k;

  tempVector (valu, bpts);
  tempVector (mult, bpts);
  tempVector (tmp , nb);

  dzero (bpts, valu, 1);
  dzero (bpts, mult, 1);

  for(k = 0, fptr = field; k < nel; k++, fptr += ntot ) {
    dgathr(nb, fptr, U[k].emap, tmp);
    for(i = 0; i < nb; i++) {
      valu[ bmap[k][i] ] += tmp[i];
      mult[ bmap[k][i] ] += 1.;
    }
  }
  
  dvdiv (bpts, valu, 1, mult, 1, valu, 1);
  
  for(k = 0, fptr = field; k < nel; k++, fptr += ntot) {
    dgathr(nb, valu, bmap[k], tmp);
    dscatr(nb, tmp , U[k].emap, fptr);
  }

  freeVector(valu); 
  freeVector(mult); 
  freeVector(tmp);
  return;
}
