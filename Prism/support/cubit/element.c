/* Element implementation
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "veclib/veclib.h"

#include "cubit.h"
#include "element.h"
#include "vertex.h"
#include "edge.h"
#include "curve.h"
#include "isomesh.h"

/* ------------------------------------------------------------------------ *
 * build_emap() - Compute the Elemental Boundary/Interior Scatter Vector    *
 *                                                                          *
 * This function computes the index map which takes a sequentially ordered  *
 * matrix into a boundary/interior ordered matrix.  The map is recomputed   *
 * whenever this function is called with a new value of (nr,ns); otherwise, *
 * the same map computed for the previous (nr,ns) is returned.              *
 * ------------------------------------------------------------------------ */

static int *build_emap (int nr, int ns)
{
  int pos  = 0;
  static int *emap = NULL, _nr = 0, _ns = 0;
  int i, j, k;

  if (emap == NULL || nr != _nr || ns != _ns) {
    _nr  = nr;
    _ns  = ns;
    emap = (int*) calloc (nr * ns, sizeof(int));
    assert (emap != NULL);

    /* Boundary scatter vector */

    for(j = 0, i = 0; j < nr; j++)
      emap[pos++] = i * nr + j;
    for(j = nr-1, i = 1; i < ns; i++)
      emap[pos++] = i * nr + j;
    for(j = nr-2, i = ns-1; j >= 0; j--)
      emap[pos++] = i * nr + j;
    for(j = 0, i = ns-2; i > 0;i--)
      emap[pos++] = i * nr + j;
    
    /* Interior scatter vector */
    
    for (i = 1;i < ns - 1; i++)
      for (j = 1;j < nr - 1; j++)
	emap[pos++] = i * nr + j;
  }
    
  return emap;
}

/* Allocate memory for a multi-frame field */

static double **frame_alloc (int nr, int ns, int nz, int nel)
{
  double *dblock, **pblock;
  int i;

  int ntot = nr * ns * nz * nel;  /* Total storage for all frames */

  /* Parallel Transpose:  Make sure the total storage area can be //
  // be partitioned (N N K) (M/P) or (N N K / P) (M).             */

  if (nz > 1) {
    const int nprocs = option("nprocs");
    ntot  = ((ntot + nprocs - 1) / nprocs) * nprocs;
  } 

  pblock = (double**) malloc (ns * nz * nel * sizeof(double*));
  dblock = (double* ) calloc (ntot, sizeof(double));
  assert (pblock != NULL && dblock != NULL);
  
  for (i = 0; i < ns * nel * nz; i++)
    pblock[i] = dblock + i * nr;

  return pblock;
}

/* ------------------------------------------------------------------------- */

Element *Element_alloc (int nr, int ns, int nz)
{
  Element *elmt = (Element*) calloc (1, sizeof(Element));
  int np = MAX(nr,ns);
  int i;

  assert(elmt);

  if (nr != ns)
    fprintf (stderr, "warning: (r,s)-polynomial order must be the same\n");
  
  /* Initialize vertices and edges */

  elmt->id        = -1;
  elmt->family    = -1;
  elmt->key       = -1;
  elmt->nr        = np;
  elmt->ns        = np;
  elmt->nz        = nz;
  elmt->emap      = build_emap  (np, np);
  elmt->xmesh     = frame_alloc (np, np, 2, 1);
  elmt->ymesh     = elmt->xmesh + np;
  elmt->vert_list = Vertex_valloc(np);
  elmt->edge_list = Edge_valloc(np, elmt->vert_list);

  return elmt;
}

void Element_free (Element *elmt)
{
  Edge_vfree  (elmt->edge_list);
  Vertex_vfree(elmt->vert_list);
  free (elmt);
}

/* ------------------------------------------------------------------------- */

void Element_setGeometry (Element *elmt, double *xc, double *yc)
{
  int i;
  
  for (i = 0; i < 4; i++) {
    ((*elmt->xmesh)[elmt->edge_list[i].start]) = xc[i];
    ((*elmt->ymesh)[elmt->edge_list[i].start]) = yc[i];
  }

  Element_genxy(elmt);
}

/* Regenerate the local mesh */

void Element_genxy (Element *elmt) 
{
  shape   (elmt);
  blend   (elmt);
  map     (elmt);
  normals (elmt);
}

void Element_dx (Element *elmt, double *u, double *dx)
{
  const int nr    = elmt->nr;
  const int ns    = elmt->ns;
  const int nrns  = nr * ns;
  double  **dr, **ds;
  int i;

  tempVector (Ur, nrns);
  tempVector (Us, nrns);
  
  getops (nr, 0, 0, &dr, 0);
  getops (ns, 0, 0, &ds, 0);

  /* Compute (r,s) partials */
  
  dgemm ('T','N', nr,ns,nr, 1., *dr, nr,   u, nr, 0., Ur, nr);
  dgemm ('N','N', nr,ns,ns, 1.,   u, nr, *ds, ns, 0., Us, nr);

  /* Compute -X- derivatives */

  if (elmt->sx)
    for (i = 0; i < nrns; i++)
      dx[i] = Ur[i] * (*elmt->rx)[i] + Us[i] * (*elmt->sx)[i];
  else
    for (i = 0; i < nrns; i++)
      dx[i] = Ur[i] * (*elmt->rx)[i];
  
  freeVector (Ur); 
  freeVector (Us);
}

void Element_dy (Element *elmt, double *u, double *dy)
{
  const int nr    = elmt->nr;
  const int ns    = elmt->ns;
  const int nrns  = nr * ns;
  double  **dr, **ds;
  int i;

  tempVector (Ur, nrns);
  tempVector (Us, nrns);
  
  getops (nr, 0, 0, &dr, 0);
  getops (ns, 0, 0, &ds, 0);

  /* Compute (r,s) partials */
  
  dgemm ('T','N', nr,ns,nr, 1., *dr, nr,   u, nr, 0., Ur, nr);
  dgemm ('N','N', nr,ns,ns, 1.,   u, nr, *ds, ns, 0., Us, nr);

  /* Compute -Y- derivatives */

  if (elmt->ry)
    for (i = 0; i < nrns; i++)
      dy[i] = Ur[i] * (*elmt->ry)[i] + Us[i] * (*elmt->sy)[i];
  else
    for (i = 0; i < nrns; i++)
      dy[i] = Us[i] * (*elmt->sy)[i];
  
  freeVector (Ur); 
  freeVector (Us);
}

/* ------------------------------------------------------------------------- *
 * Element_project(), Element_project_2d()                                   * 
 *                                                                           * 
 * These functions interpolate data from a parent element one of its four    * 
 * children.  The argument 'quad' specifies which child the interpolation    * 
 * is for.  The only difference between _project and _project_2d is that the * 
 * latter function only interpolates the 0th frame of data.                  * 
 *                                                                           * 
 * The input arrays "u" (the parent) and "v" (the child) should both be of   * 
 * size nr*ns*nz.                                                            * 
 * ------------------------------------------------------------------------- */

void Element_project (Element *elmt, const double *u, double *v, int quad) 
{
  const int nr = elmt->nr;
  const int ns = elmt->ns;
  const int nz = elmt->nz;
  double *z, *low, *high, *rr, *ss, *hr, hs;
  int i, j, k, p, q;

  quad &= 0x3;               /* Only 2 bits, so you can pass a key */
  low   = dvector (0, nr);   /* New locations */
  high  = dvector (0, nr);
  hr    = dvector (0, nr);   /* Interpolants  */

  getops (nr, &z, 0, 0, 0);

  for (i = 0; i < nr; i++) {       
    high[i] = (1. + z[i]) * .5;    /* "high" maps (0,1) */
    low [i] = high[i] - 1.;        /* "low" maps (-1,0) */
  }
  
  rr = (quad == 0 || quad == 3) ? low : high;
  ss = (quad == 0 || quad == 1) ? low : high;

  dzero (nr*ns*nz, v, 1);

  for (k = 0; k < nz; k++) {
    for (j = 0; j < ns; j++) {
      for (i = 0; i < nr; i++) {
	for (q = 0; q < nr; q++)               /* sum factorization */
	  hr[q] = hgll (q, rr[i], z, nr);      
	for (p = 0; p < ns; p++) {
	  hs    = hgll(p, ss[j], z, ns); 
	  v[i+nr*(j+ns*k)] += hs*ddot(nr, u + nr*(p+ns*k), 1, hr, 1);
	}
      }
    }
  }

  free (low);
  free (high);
  free (hr);
}

void Element_project_2d (Element *elmt, const double *u, double *v, int quad) 
{
  const int nr = elmt->nr;
  const int ns = elmt->ns;
  double *z, *low, *high, *rr, *ss, *hr, hs;
  int i, j, k, p, q;

  quad &= 0x3;               /* Only 2 bits, so you can pass a key */
  low   = dvector (0, nr);   /* New locations */
  high  = dvector (0, nr);
  hr    = dvector (0, nr);   /* Interpolants  */

  getops (nr, &z, 0, 0, 0);

  for (i = 0; i < nr; i++) {       
    high[i] = (1. + z[i]) * .5;    /* "high" maps (0,1) */
    low [i] = high[i] - 1.;        /* "low" maps (-1,0) */
  }
  
  rr = (quad == 0 || quad == 3) ? low : high;
  ss = (quad == 0 || quad == 1) ? low : high;

  dzero (nr*ns, v, 1);

  for (j = 0; j < ns; j++) {
    for (i = 0; i < nr; i++) {
      for (q = 0; q < nr; q++)               /* sum factorization */
	hr[q] = hgll (q, rr[i], z, nr);      
      for (p = 0; p < ns; p++) {
	hs    = hgll(p, ss[j], z, ns); 
	v[i+nr*j] += hs*ddot(nr, u + nr*p, 1, hr, 1);
      }
    }
  }

  free (low);
  free (high);
  free (hr);
}
