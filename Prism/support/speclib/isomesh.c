/*
 * Isoparametric Mesh Generation
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * Notes
 * -----
 * The following routines generate the collocation points for isoparametric
 * elements with straight or curved sides.  The supported types are similar
 * to those used in NEKTON but the behavior is not entirely the same.  The 
 * following list summarizes each curve type and the information that it 
 * expects to be provided in the input file:
 *
 * Type    Flag  Usage
 * --------------------------------------------------------------------
 * SinWave  W    Frequency (beta) and magnitude (alpha).  Generates a 
 *               normal sine wave along the line connecting the end-
 *               points of the edge.  For example, to specify a sine 
 *               wave with wave-number 1 and magnitude 0.5 along edge 4 
 *               of element 1, the curved side information would be:
 *
 *               4   1    1.0    0.5    W
 *
 *
 * Param.   P    Adds a parametric curve c(t) normal to the edge.  In 
 *               the string defining the curve the parameteric variable 
 *               "t" varies from 0 to 1 long the edge, and the curve is
 *               added along the outward normal.  The previous example 
 *               could have been specified as:
 *
 *               4   1  c(t) = 0.5 * sin(1.0 * PI * t)  P
 *
 *
 * Arc      C    Fit an arc with a given radius through the endpoints
 *               of the edge.  The sign of the radius determines
 *               whether the arc is convex (bows out, radius > 0) or
 *               concave (bows in, radius < 0).  To fit an arc of rad-
 *               ius 0.75 through our element's edge, enter:
 *
 *               4   1    0.75    C
 *
 *
 * File     F    Fits the element to a curve specified as (x,y)-coords.
 *               in a file.  After fitting a spline through the curve,
 *               the program finds the point where the edge would
 *               intersect the curve and adjusts the element to fit.  
 *               Any number of elements can be fit to a single curve.
 *               Example, with edge 4 of element 1 and 2 fit to the same
 *               curve:
 *
 *               4   1    mycurve.dat    F
 *               4   2    mycurve.dat    F
 *
 * --------------------------------------------------------------------
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
#include <string.h>
#include <math.h>
#include <float.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "speclib/isomesh.h"

static struct curvsort {      /* Curve generation structure */
  int        id       ;       /* Edge number                */
  CurveType  type     ;       /* Curve type                 */
  CurveInfo *info     ;       /* Pointer to curve info      */
  Edge      *edge     ;       /* Pointer to associated edge */
  Element   *elmt     ;       /* Pointer to current element */
} curvebuf [4];

/* Internal routines */

extern int    ctypecomp  (const void *p1, const void *p2);
static int    closest    (Point p, Geometry *g);
static void   blend      (Element *U),

              genStrait  (double *x, double *y, struct curvsort *p),
              genArc     (double *x, double *y, struct curvsort *p),
              genFile    (double *x, double *y, struct curvsort *p),
              genWave    (double *x, double *y, struct curvsort *p),
              genParam   (double *x, double *y, struct curvsort *p),
              bracket    (double s[], double f[], Geometry *g, Point a,
			  Vector ap),
              cvsmooth   (int nc, double *xc, double *yc, 
			          double *xp, double *yp, struct curvsort *p);

static Point  setPoint   (double x, double y),
              edgPoint   (Element *u, int id);
static Vector setVector  (Point p1, Point p2);

static double getAngle   (double s  , Geometry *g, Point a, Vector ap),
              searchGeom (Point a, Point p, Geometry *g),
              brent      (double s[], Geometry *g, Point a, Vector ap,
			  double tol);


static Geometry *lookupGeom (char *name),
                *loadGeom   (char *name);

static Geometry *geomlist;

static void geomalloc (Element *),
            geomfree  (Element *);

/* ------------------------------------------------------------------------- */

void genxy (Element *U, Curve *curv)
{
  double   *x     = U -> xmesh [0],
           *y     = U -> ymesh [0];
  register int  i;

  for (i = 0; i < 4; i++) {            /* Load up the curve buffer */
    curvebuf[i].id   = i;
    curvebuf[i].type = curv->type [i];
    curvebuf[i].info = curv->info + i;
    curvebuf[i].edge = U->edges + i;
    curvebuf[i].elmt = U;
  }

  qsort (curvebuf, 4, sizeof(struct curvsort), ctypecomp);


  /* ----------  G E N E R A T I O N ---------- */

  for (i = 0; i < 4; i++)
    switch (curvebuf[i].type) {
    case T_Strait:
      genStrait (x, y, curvebuf + i);
      break;
    case T_Arc:
      genArc    (x, y, curvebuf + i);
      break;
    case T_File:
      genFile   (x, y, curvebuf + i);
      break;
    case T_SinWave:
      genWave   (x, y, curvebuf + i);
      break;
    case T_Parametric:
      genParam  (x, y, curvebuf + i);
      break;
    default:
      speclib_error("isomesh: unknown curved-side type");
      break;
    }

  blend  (U);      /* Blend the boundary curves to get the interior mesh */

  return;
}

/* ------------------------------------------------------------------------ *
 * geomap() - Generate Isoparametric Mapping Factors                        *
 *                                                                          *
 * This function computes the geometric information for the mapping from    *
 * phsyical (x,y)-coordinates to computational (r,s)-coordinates.  The      *
 * information computed here is:                                            *
 *                                                                          *
 *          xr, yr         -  dx/dr, dy/dr                                  *
 *          xs, ys         -  dx/ds, dy/ds                                  *
 *          jac            -  xr * ys - xs * yr (Jacobian)                  *
 *          rx, ry         -  dr/dx, dr/dy                                  *
 *          sx, sy         -  ds/dx, ds/dy                                  *
 *          b              -  Mass Matrix (wr * ws * jac)                   *
 *                                                                          *
 * This is also where the family system is set up.  To disable the use of   *
 * families, make sure DisableFamilies() is called before this function!    *
 * ------------------------------------------------------------------------ */

void geomap (Element *U)
{
  const int nr    = U->nr;
  const int ns    = U->ns;
  const int nrns  = nr*ns;
  double *z, *wr, *ws, **dr, **ds;
  Edge   *edg;
  Family *fam;
  register int i, j;

  /* Check to see if this element is in a current family.  *
   * If so, just add it and return.                        */
   
  if ((fam = Family_get (U))) Family_add (fam, U);      

  /* Otherwise, we need to generate a map for this element */

  else {

    getops (nr, &z, &wr, &dr, 0);
    getops (ns, &z, &ws, &ds, 0);

    /* Now we have to do a full setup for the elemental geometry. *
     * Memory is allocated here and later reserved for other pos- *
     * sible members of the family by the call to CreateFamily(). */

    geomalloc (U);

    /* Compute the partial derivatives and Jacobian for this element */

    dgemm ('T','N', nr,ns,nr, 1., *dr,nr,*U->xmesh,nr, 0.,*U->xr,nr);
    dgemm ('N','N', nr,ns,ns, 1., *U->xmesh,nr,*ds,ns, 0.,*U->xs,nr);
    dgemm ('T','N', nr,ns,nr, 1., *dr,nr,*U->ymesh,nr, 0.,*U->yr,nr);
    dgemm ('N','N', nr,ns,ns, 1., *U->ymesh,nr,*ds,ns, 0.,*U->ys,nr);

    dvmul (nrns, *U->xr, 1, *U->ys, 1, *U->jac, 1);
    dvvtvm(nrns, *U->xs, 1, *U->yr, 1, *U->jac, 1, *U->jac, 1);
    dneg  (nrns, *U->jac,1);
 

    /* Compute the inverse partial derivatives */

    dvdiv (nrns, *U->ys, 1, *U->jac, 1, *U->rx, 1);
    dvdiv (nrns, *U->xs, 1, *U->jac, 1, *U->ry, 1);
    dvdiv (nrns, *U->yr, 1, *U->jac, 1, *U->sx, 1);
    dvdiv (nrns, *U->xr, 1, *U->jac, 1, *U->sy, 1);
    
    dneg  (nrns, *U->ry, 1);
    dneg  (nrns, *U->sx, 1);

    /* Compute the diagonal of the mass matrix */

    for (i = 0; i < ns; i++)     /* Collocate weights */
      for (j = 0; j < nr; j++)
	(U->mass)[i][j] = wr[j] * ws[i] * (U->jac)[i][j];
    
    /* Free the excess factors for non-deformed elements  *
     * (if any) and create the new family.                */

    geomfree     (U);       
    Family_create(U);
  }

  /* Compute the area array (for surface integrals) */
  
  for (edg = U->edges; edg; edg = edg->next) {
    int     np    = edg->np;
    int     start = edg->start;
    int     skip  = edg->skip;
    double *area  = edg->area;
    double *xp, *yp;
      
    if (edg->id & 1) {               /* Edge lies along r = constant */
      yp   = *U->ys + start;
      xp   = (U->xs) ? *U->xs + start : (double*) NULL;
    } 
    else {                           /* Edge lies along s = constant */
      xp   = *U->xr + start;
      yp   = (U->yr) ? *U->yr + start : (double*) NULL;
    }
      
    if (xp) dvmul (np, xp, skip, xp, skip, area, 1);
    else    dzero (np, area, 1);
    if (yp) dvvtvp(np, yp, skip, yp, skip, area, 1, area, 1);
      
    dvsqrt (np, area, 1, area, 1);
  }

  return;
}

/* ------------------------------------------------------------------------ *
 * normals() - Compute unit outward normals                                 *
 *                                                                          *
 * This function computes the unit outward normal (unx, uny) along all ele- *
 * ment boundaries.  The input is a single vector field, which is assumed   *
 * to share edge structures with all other vector fields of interest. The   *
 * computed normals are stored in "unx" and "uny" along each edge.          *
 * ------------------------------------------------------------------------ */

void normals (Element *U)
{
  int      np, start, skip, id;
  double   sign[4];
  double  *unx, *uny, len;
  Edge    *edg;
  register int i;

  const double tol = FLT_EPSILON;

  sign[0] = -1.0;         /* This is a correction to make all normals "out" */
  sign[1] =  1.0;
  sign[2] =  1.0;
  sign[3] = -1.0;
  
  for (edg = U->edges; edg; edg = edg->next) {

    id    = edg->id;
    np    = edg->np;
    unx   = edg->unx;
    uny   = edg->uny;
    start = edg->start;
    skip  = edg->skip;
    
    switch (id) {

    case 0: case 2:     /* Bottom and Top edges */

      dsmul (np, sign[id], *U->sy + start, skip, uny, 1);
      if (U->sx) 
	dsmul(np, sign[id], *U->sx + start, skip, unx, 1);
      else
	dzero(np, unx, 1);
      break;

    case 1: case 3:     /* Left and Right edges */

      dsmul(np, sign[id], *U->rx + start, skip, unx, 1);
      if (U->ry)
	dsmul(np, sign[id], *U->ry + start, skip, uny, 1);
      else
	dzero(np, uny, 1);
      break;
    }

    /* Re-scale the normals to make them "unit" normals */

    for (i = 0; i < np; i++) {
      len     = sqrt(unx[i] * unx[i] + uny[i] * uny[i]);
      unx[i] /= len;
      uny[i] /= len;
      
      if (len < tol) 
	speclib_error("Normal computation failed...check your mesh");
    }
  }

  return;
}

/*
 * Generate a Z-mesh 
 */

double *zmesh (int nz)
{
  double *z = dvector(0, nz);

  z[0] = 0.;
  z[1] = scalar("2*PI/BETA") / nz;
  dramp (nz, z+1, z+1, z+1, 1);

  return z;
}

/* Allocate memory for an Element's geometry */

static void geomalloc (Element *U)
{
  const int nrm1 = U->nr - 1;
  const int nsm1 = U->ns - 1;

  U -> rx   = dmatrix (0, nsm1, 0, nrm1);
  U -> ry   = dmatrix (0, nsm1, 0, nrm1);
  U -> sx   = dmatrix (0, nsm1, 0, nrm1);
  U -> sy   = dmatrix (0, nsm1, 0, nrm1);

  U -> jac  = dmatrix (0, nsm1, 0, nrm1);
  U -> mass = dmatrix (0, nsm1, 0, nrm1);
  
  U -> xr   = dmatrix (0, nsm1, 0, nrm1);
  U -> xs   = dmatrix (0, nsm1, 0, nrm1);
  U -> yr   = dmatrix (0, nsm1, 0, nrm1);
  U -> ys   = dmatrix (0, nsm1, 0, nrm1);

  return;
}

/* Free the excess memory for a non-deformed element */

static void geomfree (Element *U)
{
  int    nrns = U->nr * U->ns;
  double tol  = dparam("TOLABS");

  if (dasum(nrns, *U->sx, 1) < tol &&
      dasum(nrns, *U->xs, 1) < tol ) {
     free_dmatrix (U->sx, 0, 0); U->sx = (double **) NULL;
     free_dmatrix (U->xs, 0, 0); U->xs = (double **) NULL;
  }
      
  if (dasum(nrns, *U->ry, 1) < tol &&
      dasum(nrns, *U->yr, 1) < tol ) {
     free_dmatrix (U->ry, 0, 0); U->ry = (double **) NULL;
     free_dmatrix (U->yr, 0, 0); U->yr = (double **) NULL;
  }

  return;
}

/* ------------------------------------------------------------------------ */

/* Compute the remaining interior points as a bilinear blend of the  *
 * functions describing the boundaries.                              */

static void blend (Element *U)
{
  const int nr = U->nr;
  const int ns = U->ns;
  double**  x  = U->xmesh;
  double**  y  = U->ymesh;

  double   *zr, *zs;
  register int i, j;

  getops (nr, &zr, 0, 0, 0);
  getops (ns, &zs, 0, 0, 0);

  for (i = 1; i < ns-1; i++) { 
    for (j = 1; j < nr-1; j++) {
      x[i][j] = 0.50 * ((1. - zr[j])*x[i][0] + (1. + zr[j])*x[i][nr-1]
	             +  (1. - zs[i])*x[0][j] + (1. + zs[i])*x[ns-1][j])
	
	      - 0.25 * ((1. - zr[j])*(1. - zs[i])*x[0][0]
	             +  (1. - zr[j])*(1. + zs[i])*x[ns-1][0]
	             +  (1. - zs[i])*(1. + zr[j])*x[0][nr-1] 
	             +  (1. + zs[i])*(1. + zr[j])*x[ns-1][nr-1]);
    
      y[i][j] = 0.50 * ((1. - zr[j])*y[i][0] + (1. + zr[j])*y[i][nr-1]
	             +  (1. - zs[i])*y[0][j] + (1. + zs[i])*y[ns-1][j])

	      - 0.25 * ((1. - zr[j])*(1. - zs[i])*y[0][0]
	             +  (1. - zr[j])*(1. + zs[i])*y[ns-1][0]
	             +  (1. - zs[i])*(1. + zr[j])*y[0][nr-1] 
	             +  (1. + zs[i])*(1. + zr[j])*y[ns-1][nr-1]);
    }
  }

  return;
}

int ctypecomp (const void *p1, const void *p2)
{
  int t1 = (int) ((struct curvsort *) p1)->type;
  int t2 = (int) ((struct curvsort *) p2)->type;

  if (t1 < t2) return   1;
  if (t1 > t2) return  -1;
  
  return 0;
}

/* ------------------------------------------------------------------------- */

static void genArc (double *x, double *y, struct curvsort *p)
{
  Point    p1, p2, rc, rd;
  double  *xi, alpha, theta, rad, arclen;
  register int i, start, skip, np;

  genSetup (p1, p2, p, np, start, skip);
  getops   (np, &xi, 0, 0, 0);

  x     += start;
  y     += start;
  arclen = distance (p1, p2);

  if (fabs(alpha = p->info->arc.radius / arclen) < 0.505)
    speclib_error("arc radius is too small to fit an edge");
      
  alpha = alpha > 0. ? -sqrt (alpha*alpha-0.25) * arclen : 
                        sqrt (alpha*alpha-0.25) * arclen ;

  rc.x  = (p1.x + p2.x) * .5 + alpha * (p2.y - p1.y) / arclen;
  rc.y  = (p1.y + p2.y) * .5 + alpha * (p1.x - p2.x) / arclen;

  alpha = fabs (alpha);
  rad   = fabs (p->info->arc.radius);
  rd.x  = ((p1.x + p2.x) * .5 - rc.x) * (rad / alpha);
  rd.y  = ((p1.y + p2.y) * .5 - rc.y) * (rad / alpha);    
  theta = p->info->arc.radius > 0. ?  atan (.5 * arclen / alpha) :
                                     -atan (.5 * arclen / alpha) ;

  for (i = 0; i < np; i++, x += skip, y += skip) {
    alpha = theta * xi [i];
    *x    = rc.x + rd.x * cos(alpha) - rd.y * sin(alpha);
    *y    = rc.y + rd.x * sin(alpha) + rd.y * cos(alpha);
  }

  p->info->arc.xc = rc.x;
  p->info->arc.yc = rc.y;

  return;
}

static void genStrait (double *x, double *y, struct curvsort *p)
{
  Point    p1, p2;
  Vector   vedg;
  double   *xi;
  register int i, start, skip, np;

  genSetup (p1, p2, p, np, start, skip);
  getops   (np, &xi, 0, 0, 0);

  x   += start;
  y   += start;
  vedg = setVector(p1, p2);

  for (i = 0; i < np; i++, x += skip, y += skip) {
    *x = p1.x + vedg.x * .5 * (xi[i] + 1.);
    *y = p1.y + vedg.y * .5 * (xi[i] + 1.);
  }
  
  return;
}


static void genFile (double *x, double *y, struct curvsort *p)
{
  Geometry *g;
  Point    p1, p2, a;
  double   *z, *eta, xoff, yoff;
  register int i, start, skip, np;

  genSetup (p1, p2, p, np, start, skip);
  getops   (np, &z, 0, 0, 0);

  eta    = dvector (0, np);
  if ((g = lookupGeom (p->info->file.name)) == (Geometry *) NULL)
       g = loadGeom   (p->info->file.name);

  /* If the current edge has an offset, apply it now */

  xoff = p->info->file.xoffset;
  yoff = p->info->file.yoffset;
  if (xoff != 0. || yoff != 0.) {
    dsadd (g->npts, xoff, g->x, 1, g->x, 1);
    dsadd (g->npts, yoff, g->y, 1, g->y, 1);
  }

  /* Find the location along the curve where the edge intersects it */

  a         = edgPoint   (p->elmt, (p->id+3)%4);
  eta[0]    = searchGeom (a, p1, g);
  a         = edgPoint   (p->elmt, (p->id+2)%4);
  eta[np-1] = searchGeom (a, p2, g);

  /* Now generate the points where we'll evaluate the geometry */

  for (i = 1; i < np-1; i++)
    eta [i] = eta[0] + .5 * (eta[np-1] - eta[0]) * (z[i] + 1.);
  for (i = 0, x += start, y += start; i < np; i++, x += skip, y += skip) {
    *x = splint (g->npts, eta[i], g->arclen, g->x, g->sx);
    *y = splint (g->npts, eta[i], g->arclen, g->y, g->sy);
  }

  g->pos = 0;     /* reset the geometry */
  if (xoff != 0.) 
    dvsub (g->npts, g->x, 1, &xoff, 0, g->x, 1);
  if (yoff != 0.) 
    dvsub (g->npts, g->y, 1, &yoff, 0, g->y, 1);

  free (eta);    /* free the workspace */

  return;
}

static void genWave (double *x, double *y, struct curvsort *p)
{
  const int nc = 25;
  double  a  = p->info->wave.magnitude;
  double  b  = p->info->wave.wavenum;

  double  t, *z;
  Point   p1, p2, d;
  Vector  vedg;
  register int i;

  tempVector (xp, nc); tempVector (yp, nc);

  p1   = edgPoint  (p->elmt, p->id);
  p2   = edgPoint  (p->elmt,(p->id+1)%4);
  vedg = setVector (p1, p2);
  d.x  = (p2.y - p1.y) / vedg.length;
  d.y  = (p1.x - p2.x) / vedg.length;

  coef  (nc);
  getops(nc, &z, 0, 0, 0);

  for (i = 0; i < nc; i++) {
    t      =  0.5 * (z[i] + 1.);
    xp [i] =  p1.x + vedg.x * t + d.x * a * sin(M_PI * b * t);
    yp [i] =  p1.y + vedg.y * t + d.y * a * sin(M_PI * b * t);
  }

  cvsmooth (nc, xp, yp, x, y, p);

  freeVector (xp); freeVector (yp);
  return;
}

static void genParam (double *x, double *y, struct curvsort *p)
{
  const int nc = 25;
  double  *z;
  Point   p1, p2, d;
  Vector  vedg;
  register int i;

  tempVector (xp, nc);  tempVector (yp, nc); 
  tempVector (t , nc);  tempVector (c , nc);

  p1   = edgPoint  (p->elmt, p->id);
  p2   = edgPoint  (p->elmt,(p->id+1)%4);
  vedg = setVector (p1, p2);
  d.x  = (p2.y - p1.y) / vedg.length;
  d.y  = (p1.x - p2.x) / vedg.length;
  
  coef  (nc);
  getops(nc, &z, 0, 0, 0);

  for (i = 0; i < nc; i++)
    t [i] = 0.5 * (z[i] + 1.);

  vector_def ("t", p->info->map.function);
  vector_set (nc, t, c);

  for (i = 0; i < nc; i++) {
    xp [i] = p1.x + vedg.x * t[i] + d.x * c[i];
    yp [i] = p1.y + vedg.y * t[i] + d.y * c[i];
  }

  cvsmooth (nc, xp, yp, x, y, p);

  freeVector (xp); freeVector (yp); 
  freeVector (t) ; freeVector (c) ;
  return;
}

/* Interpolate based on the arclength */

static void cvsmooth (int nc, double *xp, double *yp,
		              double *x , double *y , struct curvsort *p)
{
  const int np    = p->edge->np;
  const int skip  = p->edge->skip;
  const int start = p->edge->start;

  register int  i;
  double   *z, arclen;
  tempVector (sp, nc);

  getops (np, &z, 0, 0, 0);

  x += start; 
  y += start;

  for (i = 1, *sp = 0.; i < nc; i++) 
    sp[i] = sp[i-1] + sqrt (pow(xp[i]-xp[i-1],2.)+pow(yp[i]-yp[i-1],2.));
  arclen  = sp[i-1];

  for (i = 0; i < np; i++, x += skip, y += skip) {
    *x = dpoly (nc, .5 * arclen * (z[i] + 1.), sp, xp);
    *y = dpoly (nc, .5 * arclen * (z[i] + 1.), sp, yp);
  }
  
  freeVector (sp);
  return;
}


/* ------------------------------------------------------------------------ *
 * Geometric Arithmetic                                                     *
 *                                                                          *
 * This is a small collection of routines for manipulating points and vec-  *
 * tors (see the definition of Point and Vector above).                     *
 *                                                                          *
 * ------------------------------------------------------------------------ */

static Point setPoint (double x, double y)
{
  Point p;
  p.x = x; 
  p.y = y; 
  return p;
}

static Point edgPoint (Element *U, int id)
{
  Point p;
  const int start = U->edges[id].start;

  p.x = (*U->xmesh)[start];
  p.y = (*U->ymesh)[start];

  return p;
}

static Vector setVector (Point p1, Point p2)
{
  Vector v;

  v.x      = p2.x - p1.x;
  v.y      = p2.y - p1.y;
  v.length = sqrt (v.x*v.x + v.y*v.y);

  return v;
}

/* Compute the angle between the vector ap and the vector from a to
 * a point s on the curv.  Uses the small-angle approximation */
 
static double getAngle (double s, Geometry *g, Point a, Vector ap)
{
  Point  c;
  Vector ac;

  c  = setPoint (splint(g->npts, s, g->arclen, g->x, g->sx),
                 splint(g->npts, s, g->arclen, g->y, g->sy));
  ac = setVector(a, c);
			
  return 1. - ((ap.x * ac.x + ap.y * ac.y) / (ap.length * ac.length));
}

/* Search for the named Geometry */

static Geometry *lookupGeom (char *name)
{
  Geometry *g = geomlist;

  while (g) {
    if (strcmp(name, g->name) == 0) 
      return g;
    g = g->next;
  }

  return (Geometry *) NULL;
}

/* Load a geometry file */

static Geometry *loadGeom (char *name)
{
  Geometry* g = (Geometry *) calloc (1, sizeof(Geometry));

  char      buf [BUFSIZ];
  double    tmp[_MAX_NC];
  Point     p1, p2, p3, p4;
  FILE     *fp;  
  register int i;

  if ((fp = fopen(name, "r")) == (FILE *) NULL) {
    sprintf(buf, "couldn't find the curved-side file -- %s\n", name);
    speclib_error(buf);
  }

  while (fgets (buf, BUFSIZ, fp))    /* Read past the comments */
    if (*buf != '#') break;
  
  /* Allocate space for the coordinates */

  g -> x = (double*) calloc (_MAX_NC, sizeof(double));
  g -> y = (double*) calloc (_MAX_NC, sizeof(double));

  strcpy (g->name = (char *) malloc (strlen(name)+1), name);

  /* Read the coordinates.  The first line is already in *
   * the input buffer from the comment loop above.       */
  
  i = 0;
  while (i <= _MAX_NC && sscanf (buf,"%lf%lf", g->x + i, g->y + i) == 2) {
    i++;
    if (!fgets(buf, BUFSIZ, fp)) break;
  }
  g->npts = i;

  if (i < 2 )
    speclib_error("geometry file %s doesn't have enough points\n", g->name);
  if (i > _MAX_NC)
    speclib_error("geometry file %s has too many points\n", g->name);

  /* Allocate memory for the other quantities */

  g -> sx     = (double*) calloc (g->npts, sizeof(double));
  g -> sy     = (double*) calloc (g->npts, sizeof(double));
  g -> arclen = (double*) calloc (g->npts, sizeof(double));

  /* Compute spline information for the (x,y)-coordinates.  The vector "tmp"
     is a dummy independent variable for the function x(eta), y(eta).  */

  tmp[0] = 0.; 
  tmp[1] = 1.;
  dramp  (g->npts, tmp, tmp + 1, tmp, 1);
  spline (g->npts, 1.e30, 1.e30, tmp, g->x, g->sx);
  spline (g->npts, 1.e30, 1.e30, tmp, g->y, g->sy);

  /* Compute the arclength of the curve using 4 points per segment */

  for (i = 0; i < (*g).npts-1; i++) {
    p1 = setPoint (g->x[i], g->y[i] );
    p2 = setPoint (splint (g->npts, i+.25, tmp, g->x, g->sx),
		   splint (g->npts, i+.25, tmp, g->y, g->sy));
    p3 = setPoint (splint (g->npts, i+.75, tmp, g->x, g->sx),
		   splint (g->npts, i+.75, tmp, g->y, g->sy));
    p4 = setPoint (g->x[i+1], g->y[i+1]);

    g->arclen [i+1] = g->arclen[i] + distance (p1, p2) + distance (p2, p3) +
                                     distance (p3, p4);
  }

  /* Now that we have the arclength, compute x(s), y(s) */

  spline (g->npts, 1.e30, 1.e30, g->arclen, g->x, g->sx);
  spline (g->npts, 1.e30, 1.e30, g->arclen, g->y, g->sy);

  /* add to the list of geometries */

  g->next  = geomlist;
  geomlist = g;

  fclose (fp);
  return g;
}

/* 
 * Find the point at which a line passing from the anchor point "a" 
 * through the search point "p" intersects the curve defined by "g".
 * Always searches from the last point found to the end of the curve.
 */

static double searchGeom (Point a, Point p, Geometry *g)
{
  Vector   ap;
  double   tol = dparam("TOLCURV"), s[3], f[3];
  register int ip;

  /* start the search at the closest point */

  ap   = setVector (a, p);
  s[0] = g -> arclen[ip = closest (p, g)];
  s[1] = g -> arclen[ip + 1];

  bracket (s, f, g, a, ap);
  if (fabs(f[1]) > tol) 
    brent (s, g, a, ap, tol); 

  return s[1];
}

/* ---------------  Bracketing and Searching routines  --------------- */

static int closest (Point p, Geometry *g)
{
  const double* x = g->x    + g->pos;
  const double* y = g->y    + g->pos;
  const int     n = g->npts - g->pos;

  double   len[_MAX_NC];
  register int i;
  
  for (i = 0; i < n; i++)
    len[i] = sqrt (pow(p.x - x[i],2.) + pow(p.y - y[i],2.));

  i = idmin (n, len, 1) + g->pos;
  i = MIN (i, g->npts-2);

  /* If we found the same position and it's not the very first *
   * one, start the search over at the beginning again.  The   *
   * test for i > 0 makes sure we only do the recursion once.  */

  if (i && i == g->pos) { g->pos = 0; i = closest (p, g); }

  return g->pos = i;
}

#define GOLD      1.618034
#define CGOLD     0.3819660
#define GLIMIT    100.
#define TINY      1.e-20
#define ZEPS      1.0e-10
#define ITMAX     100

#define SIGN(a,b)     ((b) > 0. ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)  (a)=(b);(b)=(c);(c)=(d);
#define SHFT2(a,b,c)   (a)=(b);(b)=(c);

#define fa f[0]
#define fb f[1]
#define fc f[2]
#define xa s[0]
#define xb s[1]
#define xc s[2]

static void bracket (double s[], double f[], Geometry *g, Point a, Vector ap)
{
  double ulim, u, r, q, fu;

  fa = getAngle (xa, g, a, ap);
  fb = getAngle (xb, g, a, ap);

  if (fb > fa) { SHFT (u, xa, xb, u); SHFT (fu, fb, fa, fu); }

  xc = xb + GOLD*(xb - xa);
  fc = getAngle (xc, g, a, ap);

  while (fb > fc) {
    r = (xb - xa) * (fb - fc);
    q = (xb - xc) * (fb - fa);
    u =  xb - ((xb - xc) * q - (xb - xa) * r) /
              (2.*SIGN(MAX(fabs(q-r),TINY),q-r));
    ulim = xb * GLIMIT * (xc - xb);

    if ((xb - u)*(u - xc) > 0.) {      /* Parabolic u is bewteen b and c */
      fu = getAngle (u, g, a, ap);
      if (fu < fc) {                    /* Got a minimum between b and c */
	SHFT2 (xa,xb, u);
	SHFT2 (fa,fb,fu);
	return;
      } else if (fu > fb) {             /* Got a minimum between a and u */
	xc = u;
	fc = fu;
	return;
      }
      u  = xc + GOLD*(xc - xb);    /* Parabolic fit was no good. Use the */
      fu = getAngle (u, g, a, ap);             /* default magnification. */

    } else if ((xc-u)*(u-ulim) > 0.) {   /* Parabolic fit is bewteen c   */
      fu = getAngle (u, g, a, ap);                         /* and ulim   */
      if (fu < fc) {
	SHFT  (xb, xc, u, xc + GOLD*(xc - xb));
	SHFT  (fb, fc, fu, getAngle(u, g, a, ap));
      }
    } else if ((u-ulim)*(ulim-xc) >= 0.) {  /* Limit parabolic u to the  */
      u   = ulim;                           /* maximum allowed value     */
      fu  = getAngle (u, g, a, ap);
    } else {                                       /* Reject parabolic u */
      u   = xc + GOLD * (xc - xb);
      fu  = getAngle (u, g, a, ap);
    }
    SHFT  (xa, xb, xc, u);      /* Eliminate the oldest point & continue */
    SHFT  (fa, fb, fc, fu);
  }
  return;
}

/* Brent's algorithm for parabolic minimization */

static double brent (double s[], Geometry *g, Point ap, Vector app, double tol)
{
  int    iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a  = MIN (xa, xc);               /* a and b must be in decending order */
  b  = MAX (xa, xc);
  d  = 1.;
  x  = w  = v  = xb;
  fw = fv = fx = getAngle (x, g, ap, app);

  for (iter = 1; iter <= ITMAX; iter++) {    /* ....... Main Loop ...... */
    xm   = 0.5*(a+b);
    tol2 = 2.0*(tol1 = tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {             /* Completion test */
      xb = x;
      return fx;
    }
    if (fabs(e) > tol1) {             /* Construct a trial parabolic fit */
      r = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v) * q-(x-w) * r;
      q = (q-r) * 2.;
      if (q > 0.) p = -p;
      q = fabs(q);
      etemp=e;
      e = d;

      /* The following conditions determine the acceptability of the    */
      /* parabolic fit.  Following we take either the golden section    */
      /* step or the parabolic step.                                    */

      if (fabs(p) >= fabs(.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d = CGOLD * (e = (x >= xm ? a-x : b-x));
      else {
	d = p / q;
	u = x + d;
	if (u-a < tol2 || b-u < tol2)
	  d = SIGN(tol1,xm-x);
      }
    } else
      d = CGOLD * (e = (x >= xm ? a-x : b-x));

    u  = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu = getAngle(u,g,ap,app);                     

    /* That was the one function evaluation per step.  Housekeeping... */

    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      SHFT(v ,w ,x ,u );
      SHFT(fv,fw,fx,fu)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v  = w;
	w  = u;
	fv = fw;
	fw = fu;
      } else if (fu <= fv || v == x || v == w) {
	v  = u;
	fv = fu;
      }
    }
  }                        /* .......... End of the Main Loop .......... */
  
  speclib_error("isomesh: too many iterations in brent()");
  xb = x;
  return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN
#undef fa
#undef fb
#undef fc
#undef xa
#undef xb
#undef xc
