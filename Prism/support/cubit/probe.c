/* 
 * PROBE Functions 
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
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "veclib/veclib.h"
#include "cubit.h"
#include "probe.h"

/* ---------------------------------------------------------------------- *
 * GetRS() -- Locate an (x,y)-point within an element                     *
 *                                                                        *
 * This function computes the inverse of the geometry map (r,s) > (x,y)   *
 * using a Newton-Rhapson iteration to find the "roots" of the function   *
 *                                                                        *
 *                             | X(r,s) - x |                             *
 *                    F(r,s) = |            | = 0                         *
 *                             | Y(r,s) - y |                             *
 *                                                                        *
 * where the iteration computes corrections by solving the linear system  *
 *                                                                        *
 *                            -[grad F] dz = F                            *
 *                                                                        *
 * Return value: 0 if the input is invalid                                *
 *               n if converged after n iterations                        *
 *              -n if diverged  after n iterations                        *
 * ---------------------------------------------------------------------- */
                             /* ............ Tolerances ............ */
#define MAXIT    15          /* maximum iterations for convergence   */
#define EPSILON  1.e-6       /* |x-xp| < EPSILON  error tolerance    */
#define DIVERGE  1000.       /* |x-xp| > DIVERGE  error tolerance    */
#define GETRSTOL 1.01        /* |r  s| < GETRSTOL inside an element? */

static GetRS (Element *elmt, double x, double y, double *r, double *s)
{
  double rr = *r;
  double ss = *s;

  int ok   = 0;
  int iter = 1;
  int i;

  int info;
  int ipiv[2];

  double eps;
  double dz[2], jac[2][2];

  double *zr, hr[_MAX_NORDER];
  double *zs, hs[_MAX_NORDER];

  if (elmt) {
    const int nr = elmt->nr;
    const int ns = elmt->ns;
    const double size = dsum(nr*ns, *elmt->mass, 1);

    getops(nr, &zr, 0, 0, 0);   /* Get GLL points */
    getops(ns, &zs, 0, 0, 0);
    
    rr = CLAMP(*r, -1., 1.);   /* Use the given values as a guess if they */
    ss = CLAMP(*s, -1., 1.);   /* are valid, otherwise clamp them.        */

    /* .... Start the iteration .... */
    
    do {
      
      /* Compute the interpolants */
      
      for (i = 0; i < nr; i++)
	hr[i] = hgll(i, rr, zr, nr);
      for (i = 0; i < ns; i++)
	hs[i] = hgll(i, ss, zs, ns);
      
      /* Compute the value of -F(r,s) and assign to dz */
      
      dz[0] = x; 
      dz[1] = y;
      for (i = 0; i < ns; i++) {
	dz[0] -= hs[i] * ddot(nr, hr, 1, (*elmt).xmesh[i], 1);
	dz[1] -= hs[i] * ddot(nr, hr, 1, (*elmt).ymesh[i], 1);
      }
      
      /* Check solution convergence and/or divergence */

      eps = dnrm2(2, dz, 1);
      if (eps < EPSILON * size)
	ok =  iter;
      else if (eps > DIVERGE * size)
	ok = -iter;

      /* Compute a new iterate */

      else {
	
	memset (jac, '\0', 4 * sizeof(double));     /* Clear */

	for (i = 0; i < ns; i++) {
	  jac[0][0] += hs[i] * ddot(nr, hr, 1, (*elmt).xr[i], 1);
	  jac[1][1] += hs[i] * ddot(nr, hr, 1, (*elmt).ys[i], 1);
	}

	if (elmt->xs) 
	  for (i = 0; i < ns; i++) 
	    jac[1][0] += hs[i] * ddot (nr, hr, 1, (*elmt).xs[i], 1);
	if (elmt->yr)
	  for (i = 0; i < ns; i++) 
	    jac[0][1] += hs[i] * ddot (nr, hr, 1, (*elmt).yr[i], 1);
      
	/* Solve the 2 x 2 system:  -[F'] dz = F */

	dgesv (2, 1, jac, 2, ipiv, dz, 2, info);
	
	/* Update the solution */

	rr += dz[0];
	ss += dz[1];
	
	/* Convergence test for iterates [always O(1)] */

	eps = dnrm2(2, dz, 1);
	if (eps < EPSILON)
	  ok =  iter;
	if (eps > DIVERGE)
	  ok = -iter;
      }
      
    } while (!ok && ++iter < MAXIT);

    /* ---------------  End of the Iteration --------------- */

  }
  
  *r = rr; 
  *s = ss;
  
  /* Only allow convergence if (r,s) lies in the standard domain */
  
  if (fabs(rr) > GETRSTOL || fabs(ss) > GETRSTOL) ok = -iter;
  
  return ok;
}
                
#undef  MAXIT   
#undef  EPSILON 
#undef  DIVERGE 
#undef  GETRSTOL

/* Find the closest point, return its index */

static int closest (int n, double zp, double *z)
{
  register int i;
  double dz [_MAX_NORDER];
  
  for (i = 0; i < n; i++)
    dz[i] = fabs(zp - z[i]);
  
  return idmin(n, dz, 1);
}

/* ---------------------------------------------------------------------- *
 * Probe_alloc() -- Allocate a Probe                                      *
 * ---------------------------------------------------------------------- */

#define RESET(r,s) ( r = s = 0. )

Probe* Probe_alloc (const Field *u, location_t type, double x, double y)
{
  Probe *probe = (Probe*) malloc(sizeof(Probe));

  probe->x    = x;
  probe->y    = y;
  probe->geom = (Field*) u;
  probe->elmt = FIELD_HEAD(u);

  switch (probe->type=type) {
  case PROBE_XP:
    probe->location.mesh.r  = 0.;
    probe->location.mesh.hr = dvector(0, probe->elmt->nr-1);
    probe->location.mesh.s  = 0.;
    probe->location.mesh.hs = dvector(0, probe->elmt->ns-1);
    break;
  case PROBE_INDEX:
    probe->location.node.ir = 0;
    probe->location.node.is = 0;
    break;
  }

  Probe_move(probe, x, y);

  return probe;
}

/* release the memory occupied by a Probe */

void Probe_free (Probe *p) 
{
  if (p->type == PROBE_XP) {
    free_dvector (p->location.mesh.hr, 0);
    free_dvector (p->location.mesh.hs, 0);
  }
  free (p);
}

/* Move the probe to a new location */

int Probe_move (Probe *p, double x, double y)
{
  const int nr = p->elmt->nr;
  const int ns = p->elmt->ns;
  int done = 0;

  double *zr, *zs, r, s;

  /* Initialization */

  getops (nr, &zr, 0, 0, 0);
  getops (ns, &zs, 0, 0, 0);

  switch (p->type) {
  case PROBE_XP:
    r = p->location.mesh.r;
    s = p->location.mesh.s;
    break;
  case PROBE_INDEX:
    r = zr[p->location.node.ir];
    s = zs[p->location.node.is];
    break;
  }

  /* Try locating the new point within the current element. If GetRS()  *
   * returns with status > 0 then we're all done.  Otherwise we need to *
   * do a full search over the mesh.                                    */

  if (GetRS(p->elmt, x, y, &r, &s) > 0)
    done = 1;

  if (!done) {
    Element *elmt = p->elmt->next;

    while (!done && elmt != NULL) {
      r = 0.;
      s = 0.;
      if (GetRS(elmt, x, y, &r, &s) > 0)
	done = 1;
      else
	elmt = elmt->next;
    }

    if (!done) {
      elmt = FIELD_HEAD(p->geom);
      while (!done && elmt != p->elmt) {
	r = 0.;
	s = 0.;
	if (GetRS(elmt, x, y, &r, &s) > 0)
	  done = 1;
	else
	  elmt = elmt->next;
      }
    }

    if (!done)
      return -1;
    else
      p->elmt = elmt;
  }

  /* Finish updating the probe */

  switch (p->type) {
  case PROBE_XP: {
    int i, j;

    p->location.mesh.r = r;
    p->location.mesh.s = s;

    for (i = 0; i < nr; i++)
      p->location.mesh.hr[i] = hgll(i, r, zr, nr);
    for (i = 0; i < ns; i++)
      p->location.mesh.hs[i] = hgll(i, s, zs, ns);
      
    for (p->x = 0., i = 0; i < ns; i++)
      for (j = 0; j < nr; j++)
	p->x += p->elmt->xmesh[i][j] *
	  ( (p->location.mesh.hs)[i] * (p->location.mesh.hr)[j] );

    for (p->y = 0., i = 0; i < ns; i++)
      for (j = 0; j < nr; j++)
	p->y += p->elmt->ymesh[i][j] *
	  ( (p->location.mesh.hs)[i] * (p->location.mesh.hr)[j] );
    break;
  }

  case PROBE_INDEX: {
    int ir, is;
      
    p->location.node.ir = (ir = closest (nr, r, zr));
    p->location.node.is = (is = closest (ns, s, zs));
      
    p->x = p->elmt->xmesh[is][ir];
    p->y = p->elmt->ymesh[is][ir];
    break;
  }

  default:
    break;
  }

  return 0;
}

/* Evaluate a Probe */

double Probe_eval (const Probe *p, const Field *u)
{
  int i;

  const int nr = FIELD_NR(u);
  const int ns = FIELD_NS(u);
  const int id = p->elmt->id;
  double val;

  switch (p->type) {
  case PROBE_XP:
    val = 0.;
    for (i = 0; i < ns; i++)
      val += (p->location.mesh.hs)[i] * 
	ddot(nr, p->location.mesh.hr, 1, 
	     FIELD_DATA(u,p->elmt->id)+i*nr, 1);
    break;

  case PROBE_INDEX: {
    const int index = nr*(p->location.node.is)+p->location.node.ir;
    double *uptr = FIELD_DATA(u,p->elmt->id);
    val = uptr[index];
    break;
  }
  }

  return val;
}
