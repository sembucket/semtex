/*
 * Residual Evaluation & Masking
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

#include "cubit.h"
#include "veclib/veclib.h"

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
 *                  tau = tol * max (norm(f), norm(R))                      *
 *                                                                          *
 * where "tol" should generally be in the range 1.e-6 to 1.e-8.             *
 * ------------------------------------------------------------------------ */

double Form_RHS 
  (Field *F, Mesh *mesh, Matrix *M, void (*A)(), double *g, double *r)
{
  const int nrns = MESH_NR(mesh)*MESH_NS(mesh);
  const int n    = M->ipts + M->bpts;

  const int *solve = M->solve;
  const int  nvert = M->nvert;

  double tau, eps;
  double buf[_MAX_NORDER*_MAX_NORDER];
  double *wr, *ws;
  int i, p;
  Element *elmt;
  Edge    *edge;

  /* Initialize the residual, r = -A(g,w) */
  
  A    (n, g, r);
  dneg (n, r, 1);

  /* Initialize quadrature weights */

  getops (MESH_NR(mesh), 0, &wr, 0, 0);
  getops (MESH_NS(mesh), 0, &ws, 0, 0);

  /* Natural (Nemann) boundary conditions */

  for (elmt = MESH_HEAD(mesh); elmt; elmt = elmt->next) {
    for (edge = elmt->edge_list; edge; edge = edge->next) {
      if (edge->bc && solve[nvert+EDGE_KEY(edge)]) {

	const int np = edge->np;
	const int n  = np-1;
	double *area = edge->area;
	double *w    = edge->id & 1 ? wr : ws;
	double *h    = buf;
	double *rptr = r+nvert+EDGE_KEY(edge)*(np-2)-1;
	Vertex *a    = edge->a;
	Vertex *b    = edge->b;

	if (user_bc) {
	  BC *bc    = edge->bc;
	  int class = (*user_bc)(bc, elmt, -1, -1, NULL);
	  if (class == BC_NATURAL) {
	    (*user_bc)(bc, elmt, a->id,       -1, h);
	    (*user_bc)(bc, elmt,    -1, edge->id, h+1);
	    (*user_bc)(bc, elmt, b->id,       -1, h+n);

	    r[a->key] += w[0] * area[0] * h[0];
	    for (i = 1; i < np-1; i++)
	      rptr[i] += w[i] * area[i] * h[i];
	    r[b->key] += w[n] * area[n] * h[n];
	  }
	}
      }
    }
  }

  /* Integrate the force and compute its approximate L2 norm */

  tau = 0.;
  for (elmt = MESH_HEAD(mesh); elmt; elmt = elmt->next) {
    const int k = ELEMENT_ID(elmt);
    for (i = 0; i < nrns; i++) {
      buf[i]    = (eps = (FIELD_DATA(F,k)[i] * (*elmt->mass)[i]));
      tau      +=  eps * (FIELD_DATA(F,k)[i]);
    }
    elmt_sum_data (elmt, r, buf);
  }
  
  Mask_RHS (mesh, M, r);

  tau = sqrt(tau);
  eps = sqrt(ddot(n, r, 1, r, 1));
  tau = dparam("TOLREL") * MAX (tau, eps);

  return tau;
}

/* ------------------------------------------------------------------------ *
 * Mask_RHS()                                                               *
 *                                                                          *
 * Assign a value of zero to each boundary node that corresponds to an ess- *
 * ential boundary condition.                                               *
 * ------------------------------------------------------------------------ */

void Mask_RHS (Mesh *mesh, Matrix *A, double *R)
{
  const int nvert = A->nvert;
  const int nedge = A->nedge;
  const int np    = mesh->nr - 2;
  int *solve = A->solve;
  int i, j;

  for (i = 0; i < nvert; i++)
    if (solve[i] == 0)
      R[i] = 0.;

  solve += nvert;
  R     += nvert;
  for (i = 0; i < nedge; i++)
    if (solve[i] == 0)
      dzero (np, R + i*np, 1);
  
  return;
}


/* ------------------------------------------------------------------------ *
 * Form_G()                                                                 *
 *                                                                          *
 * Compute a function that satisfies the appropriate boundary conditions on *
 * the edge of the computational domain.   The solution elsewhere serves    *
 * as an "initial guess" for the conjugate gradient iteration.  In this     *
 * sense, CG will correct this solution so it satisfies both the boundary   *
 * conditions and the differential equation.                                *
 * ------------------------------------------------------------------------ */

void Form_G (Field *u, Mesh *mesh, Matrix *A, double *g)
{
  Element *elmt;
  Edge    *edge;
  int     nb, ni, np, nvert, *bmap, *imap, *npts, **index;
  double  xbuf[_MAX_NB], *xb, *xelmt, *xloc;

  int p;

  xb    = g;
  xelmt = g + A->bpts;

  nb    = mesh->head->nr * 4 - 4;
  ni    = SQR(mesh->head->nr-2);
  np    = mesh->head->nr - 2;
  nvert = A->nvert;

  
  for (elmt = mesh->head; elmt; elmt = elmt->next) {
  
    npts  = A->map[elmt->id].npts;
    index = A->map[elmt->id].index;
    bmap  = elmt->emap;
    imap  = elmt->emap + nb;
    xloc  = FIELD_DATA(u,elmt->id);

    /* Initialize */
    
    for (p = 0; p < ni; p++)
      xelmt[ ni*elmt->id + p ] = xloc[ imap[p] ];
    for (p = 0; p < nb; p++)
      if (npts[p] == 1)
	xb [ index[p][0] ] = xloc[ bmap[p] ];
  }

  /* Now set BC's */

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    for (edge = elmt->edge_list; edge; edge = edge->next) {
      if (edge->bc && !( A->solve[ nvert + edge->key ] )) {
	Vertex *a = edge->a;
	Vertex *b = edge->b;

	if (user_bc) {
	  BC  *bc   = edge->bc;
	  int class = (*user_bc)(bc, elmt, -1, -1, NULL);
	  if (class == BC_ESSENTIAL) {
	    (*user_bc)(bc, elmt, a->id, -1, xb + a->key);
	    (*user_bc)(bc, elmt, b->id, -1, xb + b->key);
	    (*user_bc)(bc, elmt, -1, edge->id, xb+nvert+np*edge->key);
	  }
	}
      }
    }
  }
}


