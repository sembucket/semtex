/*
 * SOLVE - Functions for solving the Helmholtz equation
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * -------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"

/* Private functions */

static void  solve_boundary    (Field *, BSystem *, double *, double);
static void  solve_interior    (Field *, BSystem *, double *);
static void  couple            (Field *, BSystem *, double *);

/* ------------------------------------------------------------------------ *
 * Solve() - Spectral Element Helmholtz Solver                              *
 *                                                                          *
 * This function uses the computed coefficient matrices to solve the two    *
 * coupled systems:                                                         *
 *                                                                          *
 *     (1) Boundary :    A*   ub = F_bnd - A_12 A_22^ F_int                 *
 *                                                                          *
 *     (2) Interior :    A_22 ui = F_int - A_21 ub                          *
 *                                                                          *
 * where A* = (A_11 - A_12 A_22^ A_21) is the statically-condensed bound-   *
 * ary matrix.                                                              *
 *                                                                          *
 *                                                                          *
 * Input:   U  ...   initial guess                                          *
 *          F  ...   force                                                  *
 *        Ubc  ...   boundary conditions                                    *
 *          B  ...   Helmholtz system matrices                              *
 *                                                                          *
 * Output:  U  ...   solution                                               *
 * ------------------------------------------------------------------------ */

void Solve (Field *U, Field *F, Bedge *Ubc, BSystem *B)
{
  Solve_A (U, F, Ubc, B, Helmholtz);
  return;
}

/* General-purpose interface */

void Solve_A (Field *U, Field *F, Bedge *Ubc, BSystem *B, Operator A)
{
  const int bpts = B->bpts;
  const int npts = B->ipts + bpts;
  const int ntot = U->nr * U->ns * B->elements;
  double tau;

  tempVector (resid, npts);
  tempVector (v,     ntot);

  /* Apply boundary conditions to U and compute the RHS-vector.  If the *
   * residual is below the tolerance returned by Form_RHS, then return  *
   * without any action on U.  Otherwise, compute the correction to U   *
   * to make the residual zero.                                         */

  BC_set (Ubc, U, B);

  if (B->singular) dzero (ntot, *U->field, 1);
  
  tau = Form_RHS (U, F, Ubc, B, A, resid);
        Mask_RHS (U, B, resid);

  if (tau < sqrt(ddot(npts, resid, 1, resid, 1))) {
    solve_boundary  (U, B, resid, tau);
    solve_interior  (U, B, resid);

    Field_local (U, B, v, resid);
    dvadd       (ntot, v, 1, *U->field, 1, *U->field, 1);
  }

  freeVector (resid);
  freeVector (v);
  return;
}


/* -----------------  P R I V A T E    F U N C T I O N S  ----------------- */

/* 
 * Solve the boundary system (semi-direct version)
 */

static void solve_boundary (Field *U, BSystem *B, double *R, double tau)
{
  const int bdof = B->bdof;

  if (bdof) {
    double   alpha, beta, eps, epso, epsr;
    int      info;
    register int k, iter;

    double   ptmp [_MAX_NB];
    double   qtmp [_MAX_NB];

    const int  nel = B->elements;
    const int  bw  = B->bandwidth;
    int        n   = B->bdof;
    Matrix_SC* SC  = B->SC;
    double*    v   = R;

    tempVector (p, n);
    tempVector (q, n);
    tempVector (r, n);

    /* Couple the RHS of the interior system to the boundary and *
     * condense the system if it's singular (i.e., pressure).    */

    couple (U, B, R);

    if (B->singular) {
      n--; 
      R[n] = 0.;
      p[n] = 0.;
      q[n] = 0.;
      r[n] = 0.;
      v[n] = 0.;
    }

    dcopy (n, R, 1, r, 1); 
    dzero (n, v, 1);

    epsr = sqrt(ddot(n, r, 1, r, 1));
    iter = 0;

    while (epsr > tau  &&  iter++ < n) 
      {
	/* Preconditioner */
	
	dcopy (n,  r,  1, q, 1);
	if (bw) dpbtrs ('U', n, bw, 1, B->Hp, bw+1, q, n, info); 
	else    dpptrs ('U', n,     1, B->Hp,       q, n, info);
	eps = ddot (n, r, 1, q, 1);

	if (iter == 1) 
	  dcopy  (n, q, 1, p, 1);
	else {
	  beta = eps / epso;
	  dsvtvp (n, beta, p, 1, q, 1, p, 1);
	}


	/* ------------ Boundary Inner Product ------------ */
	
	dzero (n, q, 1);

	for (k = 0; k < nel; k++) {
	  const int nrows = SC[k].nrows;
	  int*  rmap      = SC[k].rowmap;

	  dgathr (nrows, p, rmap, ptmp);
	  dgathr (nrows, q, rmap, qtmp);
	  dgemv  ('N', nrows, nrows, 1., SC[k].A_11, nrows, 
		  ptmp, 1, 1., qtmp, 1);
	  dscatr (nrows, qtmp, rmap, q);
	}

	if (B->singular) q[n] = 0.;

	/* ------------ Boundary Inner Product ------------ */

	alpha = eps / ddot(n, p, 1, q, 1);
	epso  = eps;
	
	daxpy (n, -alpha, q, 1, r, 1);
	daxpy (n,  alpha, p, 1, v, 1);

	epsr  = sqrt(ddot(n, r, 1, r, 1));
      }
    
    freeVector (r);
    freeVector (q);
    freeVector (p);
  }
  
  return;
}

/*
 * compute the solution on the interior of the elements
 */

static void solve_interior (Field *U, BSystem *B, double *R)
{
  const int  nel = B->elements;
  const int  ni  = (U->nr - 2)*(U->ns - 2);
  Matrix_SC* SC  = B->SC;
  
  double* v  = R;
  double* vi = R + B->bpts;
  double  vb [_MAX_NB];

  int  k;
  for (k = 0; k < nel; k++, vi += ni) {
    const int* rmap  = SC[k].rowmap;
    const int  nrows = SC[k].nrows;
    int   i, info;
    
    for (i = 0; i < nrows; i++) 
      vb[i] = v [rmap[i]];
    
    dpptrs('U', ni, 1,          SC[k].A_22, vi, ni, info);
    dgemv ('N', ni, nrows, -1., SC[k].A_12, ni, vb, 1,  1., vi, 1);
  }

  return;
}

/*
 * boundary-interior coupling 
 */

static void couple (Field *U, BSystem *B, double *R)
{
  const int  nel = B->elements;
  const int  ni  = (U->nr - 2)*(U->ns - 2);
  Matrix_SC* SC  = B->SC;

  double* F   = R;
  double* Fi  = R + B->bpts;
  double  Fb [_MAX_NB];

  int  k;
  for (k = 0; k < nel; k++, Fi += ni) {
    const int* rmap  = SC[k].rowmap;
    const int  nrows = SC[k].nrows;
    int   i;

    for (i = 0; i < nrows; i++)
      Fb [i] = F [rmap[i]];

    dgemv ('T', ni, nrows, -1., SC[k].A_12, ni, Fi, 1, 1., Fb, 1);

    for (i = 0; i < nrows; i++) 
      F [rmap[i]] = Fb [i];
  }

  return;
}
