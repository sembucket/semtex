/*
 * Spectral Element "Operators"
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * -------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>

#include "cubit.h"
#include "veclib/veclib.h"

typedef struct sp_ops {    /* ........... Spectral Operators .......... */
  int           np    ;    /* Polynomial order                          */
  double *      z     ;    /* Gauss-Lobatto Legendre collocation points */
  double *      w     ;    /* Gauss-Lobatto Legendre weights            */
  double **     d     ;    /* Associated derivative operator            */
  double **     dt    ;    /* Transpose of the derivative matrix        */
  struct sp_ops *next ;    /* Pointer to the next one                   */
} SpOps;

SpOps *sp_list = (SpOps *) NULL;

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

  if (sp == NULL) {
    coef(np);
    sp = sp_list;
  }

  if (z)  *z  = sp->z;      /* only assign the ones requested */
  if (w)  *w  = sp->w;
  if (d)  *d  = sp->d;
  if (dt) *dt = sp->dt;

  return;
}

/* Compute a Legendre transform using data defined on a GLL grid */

#define U(i,j) (u[i*m + j])
#define A(i,j) (a[i*m + j])

void legcoef (int n, int m, const double *u, double *a)
{
  const int N = n-1;
  const int M = m-1;

  int i, j, p, q;

  double *z, *rho;
  getops (n, &z, &rho, 0, 0);

  for (i = 0; i <= N; i++) {
    const double ci = (i < N) ? (i+0.5) : .5*N;
    for (j = 0; j <= M; j++) {
      const double cj = (j < M) ? (j+0.5) : .5*M;
      A(i,j) = 0.;
      for (p = 0; p < n; p++) {
	const double P = pnleg(z[p],i);
	for (q = 0; q < m; q++) {
	  A(i,j) += rho[p] * rho[q] * U(p,q) * P * pnleg(z[q],j);
	}
      }
      A(i,j) *= (ci*cj);
    }
  } 
}

#undef U
#undef A
