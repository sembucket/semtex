/*
 * Residual Evaluation & Masking
 *
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "mortar/mortar.h"

/* ------------------------------------------------------------------------ *
 * Form_RHS_q() - Compute the RHS of the algebraic system (residual)        *
 *                                                                          *
 * This routine computes the initial residual of the discrete algebraic     *
 * system, defined as:                                                      *
 *                                                                          *
 *                        R  = -(f,w) + [h,w] - A(g,w)                      *
 *                                                                          *
 * where A is the discrete weak form of a linear operator, f is the applied *
 * force, g is a given function that satisfies the essential (Dirichlet)    *
 * boundary conditions, and h is a boundary function that satisfies the     *
 * natural (Neumann) boundary conditions.                                   *
 *                                                                          *
 * Form_RHS returns the suggested tolerance "tau" for the final solution,   *
 * given by                                                                 *
 *                                                                          *
 *                  tau = tol * max ( || f ||,  || R || )                   *
 *                                                                          *
 * where "tol" should generally be in the range 1.e-6 to 1.e-8.             *
 * ------------------------------------------------------------------------ */

double Form_RHS_q
  (Field *U, Field *F, Bedge *Ubc, BSystem *B, Operator A, double *resid)
{
  const int nrns = U->nr * U->ns;
  const int nel  = B->elements;
  const int ntot = nrns * nel;
  double  *wr, *ws, tau, eps;
  Bedge   *bc;
  int i, k, p;

  Patch *P = B->other;

  getops (U->nr, 0, &wr, 0, 0);
  getops (U->ns, 0, &ws, 0, 0);
  
  /* Compute the mortar projection of the boundary conditions */

  Project_m_s (U, P, *U->field);

  /* Initialize the residual, R = -A(g,w) */

  A    (U, B, *U->field, resid);
  dneg (ntot, resid, 1);

  /* Natural (Neumann) boundary conditions
   *
   * Now go back and add the surface integrals to the RHS for each boundary
   * with Neumann data.  This is a completely local operation, done edge by
   * edge.
   * ----------------------------------------------------------------------- */

  for (bc = Ubc; bc; bc = bc->next) {
    if (BC_getType(bc->type)==NEUMANN) {
      if (isupper(bc->type)) {
	const int np    = bc->edge->np;
	const int start = bc->edge->start;
	const int skip  = bc->edge->skip;
	double  *h      = bc->bc.value + bc->elmt->frame * np;
	double *area   = bc->edge->area;
	double *w      = bc->edge->id & 1 ? ws : wr;
	
	/* NB: "p" is the absolute index of the edge node  //
        //     "i" is the local index of the edge node.    */

	p = bc->elmt->id * nrns + start;

	for (i = 0; i < np; i++, p += skip)
	  resid [p] += w[i] * area[i] * h[i];
      }
    }
  }

  /* Integrate the force and compute its approximate L2 norm */

  tau = 0.;
  for (k = 0, p = 0; k < nel; k++)
    for (i = 0; i < nrns; i++, p++) {
      resid[p] -= (eps = (*F[k].field)[i] * (*F[k].mass)[i]);
      tau      +=  eps * (*F[k].field)[i];
    }

  tau = sqrt (tau);
  eps = sqrt (ddot(ntot, resid, 1, resid, 1));
  tau = CLAMP(tau, dparam("TOLABS"), dparam("TOLREL") * MAX(tau, eps));

  Project_s_m (U, P, resid);                 /* Project to the mortars   */
  Mask_RHS_q  (U, B, resid);                 /* Apply masks and sum      */

  return tau;
}

/* 
 * Mask residuals at boundary nodes with Dirichlet data 
 */

void Mask_RHS_q (Element *U, BSystem *B, double *u)
{
  const int bpts = B->bpts;
  const int nbcs = B->bpts - B->bdof;
  const int nel  = B->elements;
  const int nb   =(U->nr + U->ns - 2) << 1;
  const int nrns = U->nr * U->ns;

  int   **bmap = B->bmap;
  double *up;
  int i, k;
  
  tempVector (valu, bpts);
  tempVector (mult, bpts);
  
  dzero (bpts, valu, 1);     /* Clear the arrays for accumulating the */
  dzero (bpts, mult, 1);     /* edge values and the multiplicity      */

  for (k = 0, up = u; k < nel; k++, up += nrns)
    for (i = 0; i < nb; i++) {
      valu [bmap[k][i]] += up[U[k].emap[i]];
      mult [bmap[k][i]] += 1.;
    }
  
  dvdiv (bpts, valu, 1, mult, 1, valu, 1);
  dzero (nbcs, valu + B->bdof, 1);

  for (k = 0, up = u; k < nel; k++, up += nrns)
    for (i = 0; i < nb; i++)
      up [U[k].emap[i]] = valu [bmap[k][i]];

  
  freeVector (mult);
  freeVector (valu);
  return;
}
