/*
 * SOLVE_Q - Functions for solving the Helmholtz equation
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
#include "mortar/mortar.h"

/* Private functions */

static void solve_boundary (Field *, BSystem *, double *, double *, double);
static void solve_interior (Field *, BSystem *, double *, double *);
static void couple         (Field *, BSystem *, double *, double *);

static void scatr_vector   (Field *, BSystem *, double *, double *, double *);
static void gathr_vector   (Field *, BSystem *, double *, double *, double *);

/* ------------------------------------------------------------------------ *
 * Solve_q() - Spectral Element Helmholtz Solver                            *
 *                                                                          *
 * This function takes the computed coefficient matrices and solves the two *
 * coupled systems:                                                         *
 *                                                                          *
 *     (1) Boundary : (A - B inv(C) trans(B)) ub = fb - B inv(C) fi         *
 *                                                                          *
 *     (2) Interior :                       C ui = fi - trans(B) ub         *
 *                                                                          *
 * where (ub|ui) is the (boundary|interior) solution at the nodes of the    *
 * spectral element mesh.                                                   *
 *                                                                          *
 * Input:   U  ...   element vector                                         *
 *         Uf  ...   right-hand side (not overwritten)                      *
 *        Ubc  ...   array of boundary conditions                           *
 *          M  ...   boundary system matrices                               *
 *                                                                          *
 * The solution is written to U.                                            *
 * ------------------------------------------------------------------------ */

void Solve_q (Field *U, Field *F, Bedge *Ubc, BSystem *B)
{
  Solve_A_q (U, F, Ubc, B, Helmholtz);
  return;
}

void Solve_A_q (Field *U, Field *F, Bedge *Ubc, BSystem *B, Operator A)
{
  const int bpts   = B->bpts;
  const int npts   = B->ipts + bpts;
  const int ntot   = U->nr * U->ns * B->elements;
  Patch*    P      = B->other;
  double    tau;

  tempVector (rb, bpts);
  tempVector (ri, npts-bpts);
  tempVector (r , ntot);

  dzero (ntot, *U->field, 1);

  /* Apply boundary conditions to U and compute the RHS-vector.  If the *
   * residual is below the tolerance returned by Form_RHS, then return  *
   * without any action on U.  Otherwise, compute the correction to U   *
   * to make the residual zero.                                         */

  BC_set  (Ubc, U, B);

  tau = Form_RHS_q (U, F, Ubc, B, A, r);

  if (tau < sqrt(ddot(ntot,r,1,r,1))) {

    scatr_vector   (U, B, r, rb, ri);
    solve_boundary (U, B,    rb, ri, tau);
    solve_interior (U, B,    rb, ri);
    gathr_vector   (U, B, r, rb, ri);

    Project_m_s    (U, P, r);    

    dvadd (ntot, r, 1, *U->field, 1, *U->field, 1);
  }

  freeVector (r);
  freeVector (rb);
  freeVector (ri);
  return;
}

/* -----------------  P R I V A T E    F U N C T I O N S  ----------------- */

/*
 * Scatter a vector "u" into boundary and interior points, "ub" and "ui".
 * This operation SUMS contributions from coincident boundary nodes, and is
 * exactly the direct stiffness summation of nodal values.
 */

static void scatr_vector 
  (Field *U, BSystem *B, double *u, double *ub, double *ui)
{
  const int nel  = B->elements,
        nb   =(U->nr + U->ns - 2) << 1,
        ni   =(U->nr - 2)*(U->ns - 2),
        nrns = U->nr * U->ns;
  register int i, k;

  dzero (B->bpts, ub, 1);

  for (k = 0; k < nel; k++, u += nrns, ui += ni) {
    const int *bmap = B->bmap[k],
          *emap = U[k].emap;

    for (i = 0; i < nb; i++)
      ub [bmap[i]] += u [emap[i]];

    emap  += nb;
    for (i = 0; i < ni; i++)
      ui [i] = u [emap[i]];
  }

  return;
}

/* Gather a vector into element form from the boundary and interior vectors */

static void gathr_vector 
  (Field *U, BSystem *B, double *u, double *ub, double *ui)
{
  const int nel  = B->elements,
        nb   =(U->nr + U->ns - 2) << 1,
        ni   =(U->nr - 2)*(U->ns - 2),
        nrns = U->nr * U->ns;
  register int i, k;

  for (k = 0; k < nel; k++, u += nrns, ui += ni) {
    const int *bmap = B->bmap[k],
          *emap = U[k].emap;

    for (i = 0; i < nb; i++)
      u [emap[i]] = ub [bmap[i]];

    emap  += nb;
    for (i = 0; i < ni; i++)
      u [emap[i]] = ui [i];
  }

  return;
}


/* 
 * Solve the boundary system (iterative version)
 */

static void solve_boundary 
   (Field *U, BSystem *B, double *ub, double *fi, double tau)
{
  int bdof = B->bdof;

  if (bdof) {
    double   alpha, beta, eps, epso, epsr;
    int      info;
    int k, iter;

    double   ptmp [_MAX_NB];
    double   qtmp [_MAX_NB];

    const int  nel = B->elements;
    const int  bw  = B->bandwidth;
    int        n   = B->bdof;
    Matrix_SC* SC  = B->SC;

    tempVector (p, n);
    tempVector (q, n);
    tempVector (r, n);
    tempVector (v, n);

    /* Couple the RHS of the interior system to the boundary and *
     * condense the system if it's singular (i.e., pressure).    */

    couple (U, B, ub, fi);

    if (B->singular) {
      --n;
      p [n] = 0.;
      q [n] = 0.;
      r [n] = 0.;
      v [n] = ub[n] = 0.;
    }

    dcopy (n, ub, 1, r, 1); 
    dzero (n, v,  1);

    epsr = sqrt(ddot(n, r, 1, r, 1));
    iter = 0;

    while (epsr > tau && iter++ < n) 
      {
	/* Preconditioner */
	
	dcopy (n, r, 1, q, 1);
#if 1
	if (bw) dpbtrs ('U', n, bw, 1, B->Hp, bw+1, q, n, info); 
	else    dpptrs ('U', n,     1, B->Hp,       q, n, info);
#else
	/* nothing: solve the identity q = q (no preconditioner) */
#endif

	eps = ddot(n, r, 1, q, 1);

	if (iter == 1) 
	  dcopy  (n, q, 1, p, 1);
	else {
	  beta = eps / epso;
	  dsvtvp (n, beta, p, 1, q, 1, p, 1);
	}


	/* ------------ Boundary Inner Product ------------ */
	
	dzero(n, q, 1);
	for (k = 0; k < nel; k++) {
	  const int nrows = SC[k].nrows;
	  int*  rmap  = SC[k].rowmap;

	  dgathr (nrows, p, rmap, ptmp);
	  dgathr (nrows, q, rmap, qtmp);

	  dgemv  ('N', nrows, nrows, 1., SC[k].A_11, nrows, ptmp, 
		  1, 1., qtmp, 1);
	  dscatr (nrows, qtmp, rmap, q);

	  if (B->singular) q[n] = 0.;
	}

	/* ------------------------------------------------ */

	alpha = eps / ddot(n, p, 1, q, 1);
	epso  = eps;
	
	daxpy (n, -alpha, q, 1, r, 1);
	daxpy (n,  alpha, p, 1, v, 1);

	epsr = sqrt (ddot(n, r, 1, r, 1));
      }

    dcopy (n, v, 1, ub, 1);
    
    freeVector (r);
    freeVector (p);
    freeVector (q);
    freeVector (v);
  }
  
  return;
}

/*
 * compute the solution on the interior of the elements
 */

static void solve_interior (Field *U, BSystem *B, double *v, double *vi)
{
  const int  nel = B->elements;
  const int  ni  = (U->nr - 2)*(U->ns - 2);
  Matrix_SC* SC  = B->SC;
  
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

static void couple (Field *U, BSystem *B, double *Fb, double *Fi)
{
  const int nel  = B->elements;
  const int ni   = (U->nr - 2)*(U->ns - 2);
  int i, k;

  Matrix_SC *SC = B->SC;

  for (k = 0; k < nel; k++, Fi += ni) {
    const int* rmap  = SC[k].rowmap;
    const int  nrows = SC[k].nrows;

    double fb [_MAX_NB];

    for (i = 0; i < nrows; i++)
      fb [i] = Fb [rmap[i]];

    dgemv ('T', ni, nrows, -1., SC[k].A_12, ni, Fi, 1, 1., fb, 1);

    for (i = 0; i < nrows; i++) 
      Fb [rmap[i]] = fb [i];
  }

  return;
}
