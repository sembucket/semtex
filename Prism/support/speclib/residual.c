/*
 * Residual Evaluation & Masking
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
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"

/* ------------------------------------------------------------------------ *
 * Form_RHS() - Compute the RHS of the algebraic system (residual)          *
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

double Form_RHS 
  (Field *U, Field *F, Bedge *Ubc, BSystem *B, Operator A, double  *R)
{
  const int nrns     = U->nr * U->ns;
  const int nel      = B->elements;
  const int ntot     = nrns * nel;
  double  *wr, *ws, tau, eps;
  Bedge   *bc;
  register int i, k, p;

  tempVector (resid, ntot);

  getops (U->nr, 0, &wr, 0, 0);
  getops (U->ns, 0, &ws, 0, 0);

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
	double  *area   = bc->edge->area;
	double  *w      = bc->edge->id & 1 ? ws : wr;
	
	/* NB: "p" is the absolute index of the edge node  
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

  tau = sqrt(tau);
  eps = sqrt(ddot(ntot, resid, 1, resid, 1));
  tau = CLAMP(tau, dparam("TOLABS"), dparam("TOLREL") * MAX(tau, eps));

  Field_dsum (U, B, resid, R);

  freeVector (resid);
  return tau;
}

/* ------------------------------------------------------------------------ *
 * Mask_RHS()                                                               *
 *                                                                          *
 * Assign a value of zero to each boundary node that corresponds to an ess- *
 * ential boundary condition.                                               *
 * ------------------------------------------------------------------------ */

void Mask_RHS (Field *U, BSystem *M, double *R)
{
  const int bpts   = M->bpts;
  const int bdof   = M->bdof;

  dzero (bpts-bdof, R + bdof, 1);

  return;
}
