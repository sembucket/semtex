/*
 * Spectral Element "Operators"
 *
 * $Id$
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 * ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "speclib/operators.h"

typedef struct sp_ops {    /* ........... Spectral Operators .......... */
  int           np    ;    /* Polynomial order                          */
  double*       z     ;    /* Gauss-Lobatto Legendre collocation points */
  double*       w     ;    /* Gauss-Lobatto Legendre weights            */
  double**      d     ;    /* Associated derivative operator            */
  double**      dt    ;    /* Transpose of the derivative matrix        */
  struct sp_ops *next ;    /* Pointer to the next one                   */
} SpOps;

SpOps* sp_list = (SpOps*) NULL;

/* ------------------------------------------------------------------------ *
 *  coef() - Generate Spectral Operators and Coefficients                   *
 *                                                                          *
 *  This function generates...                                              *
 *                                                                          *
 *                     - Derivative Operators                               *
 *                     - Interpolation Operators                            *
 *                     - Collocation Points                                 *
 *                     - Weights                                            *
 *                                                                          *
 *  ...associated with the Gauss-Lobatto Legendre mesh.                     *
 *                                                                          *
 * ------------------------------------------------------------------------ */

void coef (int norder)
{
  SpOps *sp;

  for (sp = sp_list; sp; sp = sp->next)   /* Check the current list */
    if (sp->np == norder)
      return;

  sp          = (SpOps*) malloc(sizeof(SpOps));
  sp->np      = norder;
  sp->z       = dvector(0, norder - 1);
  sp->w       = dvector(0, norder - 1);
  sp->d       = dmatrix(0, norder - 1, 0, norder - 1);
  sp->dt      = dmatrix(0, norder - 1, 0, norder - 1);
  sp->next    = sp_list;
  sp_list     = sp;

  zwgll(sp->z, sp->w, norder);
  dgll (sp->d, sp->dt, sp->z, norder);

  return;
}

void getops (int np, double *z[], double *w[], double ***d, double ***dt)
{
  SpOps *sp;

  for (sp = sp_list; sp; sp = sp->next)
    if (sp->np == np)
      break;
  if (sp == (SpOps*) NULL)
    speclib_error("invalid N requested in getops()");

  if (z)  *z  = sp->z;      /* only assign the ones requested */
  if (w)  *w  = sp->w;
  if (d)  *d  = sp->d;
  if (dt) *dt = sp->dt;

  return;
}

/* ------------------------------------------------------------------------ *
 * geofac - Geometric Factors for A u = F                                   *
 *                                                                          *
 * Compute the geometric factors (G1,G2,G3) used to evaluate the inner      *
 * product a(u,w).                                                          *
 *                                                                          *
 *          G1     =>         (xs^2 + ys^2) * wr * ws / jac                 *
 *          G2     =>         (xr^2 + yr^2) * wr * ws / jac                 *
 *          G3     =>  -(xr * sx + yr * ys) * wr * ws / jac                 *
 *                                                                          *
 * ------------------------------------------------------------------------ */

void geofac (Element *U, Matrix_IP *G)
{
  const int nr   = U->nr;
  const int ns   = U->ns;
  const int nrns = nr * ns;
  register int i;

  G -> g1 = dmatrix (0, ns-1, 0, nr-1);       /* Allocate Memory */
  G -> g2 = dmatrix (0, ns-1, 0, nr-1);
  G -> g3 = dmatrix (0, ns-1, 0, nr-1);

  if (U->xs)
    for (i = 0; i < nrns; i++) {
      (*G->g1)[i]  = (*U->xs)[i] * (*U->xs)[i] + (*U->ys)[i] * (*U->ys)[i];
      (*G->g3)[i]  = (*U->xr)[i] * (*U->xs)[i];
    }
  else {
    dvmul (nrns, *U->ys, 1, *U->ys, 1, *G->g1, 1);
    dzero (nrns, *G->g3, 1);
  }

  if (U->yr)
    for (i = 0; i < nrns; i++) {
      (*G->g2)[i]  = (*U->xr)[i] * (*U->xr)[i] + (*U->yr)[i] * (*U->yr)[i];
      (*G->g3)[i] += (*U->yr)[i] * (*U->ys)[i];
    }
  else
    dvmul (nrns, *U->xr, 1, *U->xr, 1, *G->g2, 1);

  for (i = 0; i < nrns; i++) {
    const double scale = (*U->mass)[i] / ( (*U->jac)[i] * (*U->jac)[i] );

    (*G->g1)[i] *=  scale;
    (*G->g2)[i] *=  scale;
    (*G->g3)[i] *= -scale;
  }

  if (!(U->xs || U->yr))
    { free_dmatrix (G->g3, 0, 0); G->g3 = NULL; }

  return;
}

/* Same as Field_grad(), but for a single element */

void Element_grad (const Element *U, Element *X, Element *Y)
{
  const int nr   = U->nr;
  const int ns   = U->ns;
  const int nrns = nr * ns;
  int i;

  tempVector (Ur, nrns);
  tempVector (Us, nrns);

  double  **dr, **ds;
  getops (nr, 0, 0, &dr, 0);
  getops (ns, 0, 0, &ds, 0);

  /* Compute (r,s) partials */

  dgemm ('T','N', nr,ns,nr, 1., *dr      ,nr, *U->field,nr, 0., Ur, nr);
  dgemm ('N','N', nr,ns,ns, 1., *U->field,nr, *ds      ,ns, 0., Us, nr);

  /* Compute -X- derivatives */

  if (X)
    if (X->sx)
      for (i = 0; i < nrns; i++)
	(*X->field)[i] = Ur[i] * (*X->rx)[i] + Us[i] * (*X->sx)[i];
    else
      for (i = 0; i < nrns; i++)
	(*X->field)[i] = Ur[i] * (*X->rx)[i];
  
  
  /* Compute -Y- derivatives */
  
  if (Y)
    if (Y->ry)
      for (i = 0; i < nrns; i++)
	(*Y->field)[i] = Ur[i] * (*Y->ry)[i] + Us[i] * (*Y->sy)[i];
    else
      for (i = 0; i < nrns; i++)
	(*Y->field)[i] = Us[i] * (*Y->sy)[i];
  
  freeVector (Ur); 
  freeVector (Us);
}

void Element_dx (const Element *U, Element *dUdx) 
{ Element_grad (U, dUdx, 0); }

void Element_dy (const Element *U, Element *dUdy) 
{ Element_grad (U, 0, dUdy); }
     
/* ------------------------------------------------------------------------- *
 * Spectral Element Operators                                                *
 *                                                                           *
 * The following routines are among the most important in the library, as    *
 * they define the spectral element operators: the matrices which act on the *
 * elemental data to generate the discrete Helmholtz equations.  There is a  *
 * function for the conjugate gradient solver and a function used to form    *
 * stiffness matrices for the static-condensation solver.  You can use a     *
 * different operator (variable coefficient, polar coordinates, etc.) by     *
 * providing a new operator function to the solvers.                         *
 *                                                                           *
 * Arguments:                                                                *
 *                                                                           *
 * Operator   ( Field   *U,    Template element array                        *
 *              BSystem *B,    Boundary system structure                     *
 *              double  *in,   Input array, default to U if (NULL)           *
 *              double  *out ) Output array, possible the same as "in"       *
 *                                                                           *
 * OperatorSC ( Element *U,    Template element structure (single!)          *
 *              BSystem *B,    Boundary system structure                     *
 *              double  **A,   Boundary-boundary coupling matrix             *
 *              double  **B,   Boundary-interior coupling matrix             *
 *              double  **C )  Interior-interior coupling matrix             *
 * ------------------------------------------------------------------------- */

void Helmholtz (Field *U, BSystem *B, double *in, double *out)
{
  const int nr     = U->nr;
  const int ns     = U->ns;
  const int nrns   = nr * ns;
  const int nel    = Field_count(U);
  const int ntot   = nrns * nel;
  register int i, k;

  double lambda = B->constant;
  Matrix_IP  *G = B->CG;

  double **dr, **ds, *Pr, *Ps;

  tempVector (tmp, nrns + ntot + ntot);  /* Get temporary workspace */

  Pr = tmp + nrns;
  Ps = tmp + nrns + ntot;

  getops (nr, 0, 0, &dr, 0);             /* Get derivative matrices */
  getops (ns, 0, 0, &ds, 0);   

  if (!in) in = *U->field;
  
  /* Compute (r,s)-derivatives */

  dgemm ('T','N', nr, ns*nel, nr, 1., *dr, nr, in, nr, 0., Pr, nr); 
  for (k = 0; k < nel; k++, in += nrns, Ps += nrns) 
    dgemm ('N','N', nr, ns, ns, 1., in, nr, *ds, ns, 0., Ps, nr); 

  in -= ntot;
  Ps -= ntot;

  for (k = 0; k < nel; k++, in += nrns, out += nrns, Pr += nrns, Ps += nrns)
    if (G[k].g3) {
      for (i = 0; i < nrns; i++) {
	tmp[i]  = (*G[k].g1)  [i] * Pr[i] + (*G[k].g3)[i] * Ps[i];
	Ps [i]  = (*G[k].g2)  [i] * Ps[i] + (*G[k].g3)[i] * Pr[i];
	out[i]  = (*U[k].mass)[i] * in[i];
	Pr [i]  = tmp [i];
      }
    } else {
      for (i = 0; i < nrns; i++) {
	Pr [i] *= (*G[k].g1)  [i];
	Ps [i] *= (*G[k].g2)  [i];
	out[i]  = (*U[k].mass)[i] * in[i];
      }
    }

  /* Compute final derivatives */

  out -= ntot;
  Pr  -= ntot;
  Ps  -= ntot;

  dgemm ('N', 'N', nr, ns*nel, nr, 1., *dr, nr, Pr, nr, lambda, out, nr);
  for   (k = 0; k < nel; k++, out += nrns, Ps += nrns)
    dgemm ('N', 'T', nr, ns, ns, 1., Ps, nr, *ds, ns, 1., out, nr);

  freeVector (tmp); 
  return;
}

/* ------------------------------------------------------------------------- *
 * Static-condensation form of the Helmholtz Matrix                          *
 *                                                                           *
 * NOTE: The output arrays must be dimensioned as:                           *
 *                                                                           *
 *                             A_bb  [ nb ][ nb ]                            *
 *                             A_bi  [ nb ][ ni ]                            *
 *                             A_ii  [ ni ][ ni ]                            *
 *                                                                           *
 * Also, the output matrices can be given as NULL if that section isn't      *
 * needed (i.e., when using the Family system).                              *
 * ------------------------------------------------------------------------- */

void HelmholtzSC 
  (Element *U, BSystem *B, double **A_bb, double **A_bi, double **A_ii)
{
  const double lambda = B->constant;

  const int id    = U->id;
  const int nr    = U->nr;
  const int ns    = U->ns;
  const int nrns  = nr * ns;
  const int nb    = (nr + ns - 2) << 1;
  const int ni    = (nr - 2)*(ns - 2);

  int     *emap   = U->emap;
  double  **g1    = B->CG[id].g1;
  double  **g2    = B->CG[id].g2;
  double  **g3    = B->CG[id].g3;
  double  **mass  = U->mass;
  double  **hk    = dmatrix (0, ns-1, 0, nr-1);

  double  tmp  [_MAX_NORDER*_MAX_NORDER];
  int     pmap [_MAX_NORDER*_MAX_NORDER];
  int i, j, ij, m, n;

  double  **dr, **drt, **ds, **dst;
  getops (nr, 0, 0, &dr, &drt);           /* Get derivative matrices */
  getops (ns, 0, 0, &ds, &dst);

  ij = 0;                                 /* Invert the emap array */
  for (i = 0; i < nrns; i++) 
    pmap[U->emap[i]] = ij++;

  /* .......... Compute the elemental Helmholtz matrix .......... */

    for (i = 0; i < ns; i++) {
      for (j = 0; j < nr; j++) {
	
	if ((ij = pmap[i*nr+j]) >= nb && !A_ii) 
	  continue;
	if ((ij < nb) && !(A_bb || A_bi))
	  continue;

	dzero (nrns, *hk, 1);
	
	for (m = 0; m < nr; m++) {                           
	  dvmul(nr, drt[j], 1, drt[m], 1, tmp, 1);
	  hk[i][m]  = ddot(nr, g1[i] , 1, tmp, 1);
	}
	
	for (n = 0; n < ns; n++) {                           
	  dvmul(ns, dst[i], 1, dst[n], 1, tmp, 1);
	  hk[n][j] += ddot(ns, *g2+j, ns, tmp, 1);
	}
	
	if (g3) 
	  for (n = 0; n < ns; n++) {
	    for(m = 0; m < nr; m++) {
	      hk[n][m] += g3[i][m] * dr[m][j] * ds[i][n];     
	      hk[n][m] += g3[n][j] * dr[j][m] * ds[n][i];     
	    }
	  }
	
	hk[i][j] += mass[i][j] * lambda;

	/* Copy the correct pieces to the matrices */

	if (ij < nb) {
	  if (A_bb) dgathr (nb, *hk, emap,    A_bb[ij]);     /* ....boundary */
	  if (A_bi) dgathr (ni, *hk, emap+nb, A_bi[ij]);     /* ....coupling */
	} else {
	            dgathr (ni, *hk, emap+nb, A_ii[ij-nb]);  /* ....interior */
	}
      }
    }

  free_dmatrix (hk, 0, 0);
  return;
}



