/*
 * Isoparametric Mesh Generation
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * Notes
 * -----
 * The following routines generate the mesh of internal and boundary points
 * for an element.  Internal points are blended from the edges of the 
 * element.  Edges can be one of three types: a line, a circular arc, or
 * a Bezier spline.  These are specified in the "Curves" section of the
 * input file.
 *
 * By default, the edge of an element is a straight line formed by 
 * connecting the two vertices along that edge.  For the other two types,
 * you specify the element, edge, type of curve, and additional information
 * in the "Curves" section of the mesh file.  Each line of the Curves section
 * describes one edge as follows:
 *
 *   [edge] [element] [options] [type]
 *
 * where [edge] and [element] are ID numbers, [options] are described below,
 * and [type] is L for a line, C for a circular arc, and S for a spline.
 * The syntax for each type is described next, along with an example of
 * the specification for edge 4 of element 1.
 *
 * Type    Flag  Usage
 * --------------------------------------------------------------------
 * Line     L    Default.  You don't have to specify this one, but if you
 *               really want to you just specify the type without any 
 *               additional info:
 *
 *               4  1  L
 *
 *
 * Arc      C    Fit an arc with a given radius through the endpoints of
 *               the edges.  If the radius is positive the edge will deform
 *               outwards from the center of the element, and if the 
 *               radius is negative it will deform inwards.  An arc with
 *               radius 1.2 would be given as:
 *
 *               4  1  1.2  C
 *
 *
 * Spline   S    Fit a Bezier cubic spline to the edge, defined by four
 *               control points.  Two of the control points are the endpoints
 *               of the edge, (x0,y0) and (x3,y3).  The other two are
 *               specified as [options].  The curve will be tangent to
 *               (x1,y1)-(x0,y0) at (x0,y0) and (x3,y3)-(x2,y2) at (x3,y3),
 *               and is guaranteed to stay inside the convex hull of the
 *               four points.   The spline is specified by giving the two 
 *               internal control points:
 *
 *               4  1  x1 y1 x2 y2 S
 *
 *               For example, a spline with control points (x1,y1) = (1.,1.)
 *               and (x2,y2) = (1.,0.) would be given as:
 *
 *               4  1  1.  1.  1.  0.  S
 *
 * More types may be added later.
 *
 *
 * Callable functions:
 *
 *   curve   (Element*, Edge*)        Generate a curve 
 *
 *   blend   (Element*)               Generate interior points
 *
 *   map     (Element*)               Generate an isoparametric mapping
 *
 *   normals (Element*)               Generate unit normals
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

#include "cubit.h"
#include "isomesh.h"
#include "curve.h"
#include "veclib/veclib.h"

/* ------------------------------------------------------------------------- */

/* Some functions for working with Points and Vectors in 2D */

typedef struct {
  double x;
  double y;
} Point;

typedef struct {
  double x;
  double y;
} Vector;

Point  Point_new (double x, double y) 
{ Point p; p.x = x; p.y = y; return p; }

double Point_distance (Point p1, Point p2) 
{ return sqrt(SQR(p2.x-p1.x) + SQR(p2.y-p1.y)); }
  
Point  Point_edge (Element *elmt, int id)
{ Point p; 
  p.x = (*elmt->xmesh)[elmt->edge_list[id].start];
  p.y = (*elmt->ymesh)[elmt->edge_list[id].start];
  return p; }

void Vector_setValue (Vector *v, Point p1, Point p2) { 
  v->x = p2.x-p1.x; 
  v->y = p2.y-p1.y; 
}

double Vector_length (Vector v) { 
  return sqrt(v.x*v.x + v.y*v.y); 
}

/* ------------------------------------------------------------------------- */

/* Private functions */

static void lineto   (int np, double *x, double *y, struct line   line);
static void arcto    (int np, double *x, double *y, struct arc    arc);

static void spline_fit     (spline_t *s, int np, double *x, double *y);
static void spline_compute (spline_t *s, Point, Point, Point, Point);
static void spline_getValue(spline_t *s, double t, double *x, double *y);

static double geofile_search (geofile_t *g, Point a, Point p);
static int    geofile_closest(geofile_t *g, Point p);
static void   geofile_bracket(geofile_t *g, double s[], double f[],
			      Point a, Vector ap);
static double geofile_brent  (geofile_t *g, double s[], 
			      Point a, Vector ap, double tol);
static double geofile_angle  (geofile_t *g, double s, Point a, Vector ap);
static void   geofile_fit    (geofile_t *g, Element *elmt, Edge *edge,
			      int np, double *x, double *y);
static int    geofile_reset  (geofile_t *g);

static void geomalloc (Element *);
static void geomfree  (Element *);

/* ------------------------------------------------------------------------- */

void shape (Element *elmt) 
{
  Edge *edge; 
  for(edge = elmt->edge_list; edge; edge = edge->next)
    curve (elmt, edge);
}

void curve (Element *elmt, Edge *edge)
{
  int    np    = edge->np;
  Curve *curve = edge->curve;
  double xp[_MAX_NB], yp[_MAX_NB];

  switch (curve != NULL ? curve->type : Line) {
  case Spline:
    spline_fit (&curve->info.spline, np, xp, yp);
    break;

  case Arc:
    arcto  (np, xp, yp, curve->info.arc);
    break;

  case File:
    geofile_fit (&curve->info.file, elmt, edge, np, xp, yp);
    break;

  case Line: {
    struct line line = make_line(elmt,edge);
    lineto (np, xp, yp, line);
    break;
  }

  default:
    cubit_err ("unknown curved-side type");
    break;
  }
  
  ecopy (np, xp, 1, *elmt->xmesh + edge->start, edge->skip);
  ecopy (np, yp, 1, *elmt->ymesh + edge->start, edge->skip);

  return;
}

/* ------------------------------------------------------------------------- */

void blend (Element *U)
{
  const int nr = U->nr;
  const int ns = U->ns;
  double **x = U->xmesh;
  double **y = U->ymesh;
  int i, j;

  double *zr, *zs;
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

/* ------------------------------------------------------------------------ *
 * map() - Generate Isoparametric Mapping Factors                           *
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
 * families, make sure to call Family_disable() before this function!       *
 * ------------------------------------------------------------------------ */

void map (Element *U)
{
  const int nr   = U->nr;
  const int ns   = U->ns;
  const int nrns = nr*ns;
  double *z, *wr, *ws, **dr, **ds;
  Edge   *edge;
  Family *fam;
  int i, j;

  /* Check to see if this element is in a current family.  *
   * If so, just add it and return.                        */

#if 0   
  if (fam = Family_get(U)) Family_add (fam, U);      
#else
  if (0) {
    /* nothing */
  }
#endif

  /* Otherwise, we need to generate a map for this element */

  else {

    getops (nr, &z, &wr, &dr, 0);
    getops (ns, &z, &ws, &ds, 0);

    /* Now we have to do a full setup for the elemental geometry. *
     * Memory is allocated here and later reserved for other pos- *
     * sible members of the family by the call to Family_create() */

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

    /* Compute the (diagonal) mass matrix */

    for (i = 0; i < ns; i++)     /* Collocate weights */
      for (j = 0; j < nr; j++)
	(U->mass)[i][j] = wr[j] * ws[i] * (U->jac)[i][j];
    
    /* Free the excess factors for non-deformed elements  *
     * (if any) and create the new family.                */

    geomfree     (U);       
    Family_create(U);
  }

  /* Compute the area array (for surface integrals) */
  
  for (i = 0, edge = U->edge_list; i < 4; i++) {
    int     np    = edge[i].np;
    int     start = edge[i].start;
    int     skip  = edge[i].skip;
    double *area  = edge[i].area;
    double *xp, *yp;
      
    if (edge[i].id & 1) {            /* Edge lies along r = constant */
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

static void edge_vert_init (Element *elmt)
{
  Edge   *edge  = elmt->edge_list;
  Vertex *vert  = elmt->vert_list;
  const int nb  =(elmt->nr + elmt->ns - 2) << 1;
  const int np  = edge->np;
  double pos[2], xb[_MAX_NB], yb[_MAX_NB], hh[_MAX_NORDER], *z;
  int i;

  /* We need to interpolate the midpoint coordinate */

  getops (np, &z, 0, 0, 0);
  for (i = 0; i < np; i++)
    hh[i] = hgll (i, 0., z, np);

  dgathr(nb, *elmt->xmesh, elmt->emap, xb); xb[nb] = xb[0];
  dgathr(nb, *elmt->ymesh, elmt->emap, yb); yb[nb] = yb[0];

  for (i = 0; i < 4; i++) {
    const int offset = edge[i].offset;

    pos[0] = xb[offset];
    pos[1] = yb[offset];
    memcpy (vert[i].pos, pos, 2*sizeof(double));
      
    pos[0] = ddot (np, xb + offset, 1, hh, 1);
    pos[1] = ddot (np, yb + offset, 1, hh, 1);
    memcpy (edge[i].pos, pos, 2*sizeof(double));
  }
}

void normals (Element *elmt)
{
  int np, start, skip, id;
  double sign[4];
  double *unx, *uny, len;
  Edge *edge;
  int i, j;

  const double tol = FLT_EPSILON;

  sign[0] = -1.;         /* This is a correction to make all normals "out" */
  sign[1] =  1.;
  sign[2] =  1.;
  sign[3] = -1.;
  
  edge_vert_init(elmt);  /* Initialize edge midpoints and vertices */

  for (i = 0, edge = elmt->edge_list; i < 4; i++) {

    id    = edge[i].id;
    np    = edge[i].np;
    unx   = edge[i].unx;
    uny   = edge[i].uny;
    start = edge[i].start;
    skip  = edge[i].skip;
    
    switch (id) {

    case 0: case 2:     /* Bottom and Top edges */

      dsmul (np, sign[id], *elmt->sy + start, skip, uny, 1);
      if (elmt->sx) 
	dsmul(np, sign[id], *elmt->sx + start, skip, unx, 1);
      else
	dzero(np, unx, 1);
      break;

    case 1: case 3:     /* Left and Right edges */

      dsmul(np, sign[id], *elmt->rx + start, skip, unx, 1);
      if (elmt->ry)
	dsmul(np, sign[id], *elmt->ry + start, skip, uny, 1);
      else
	dzero(np, uny, 1);
      break;
    }

    /* Re-scale the normals to make them "unit" normals */

    for (j = 0; j < np; j++) {
      len     = sqrt(unx[j] * unx[j] + uny[j] * uny[j]);
      unx[j] /= len;
      uny[j] /= len;
      
      if (len < tol) 
	cubit_err ("Normal computation failed...check your mesh");
    }
  }
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
  int nrm1 = U->nr - 1,
      nsm1 = U->ns - 1;
  
  U -> xr   = dmatrix (0, nsm1, 0, nrm1);
  U -> rx   = dmatrix (0, nsm1, 0, nrm1);
  U -> ys   = dmatrix (0, nsm1, 0, nrm1);
  U -> sy   = dmatrix (0, nsm1, 0, nrm1);

  U -> sx   = dmatrix (0, nsm1, 0, nrm1);
  U -> xs   = dmatrix (0, nsm1, 0, nrm1);
  U -> yr   = dmatrix (0, nsm1, 0, nrm1);
  U -> ry   = dmatrix (0, nsm1, 0, nrm1);

  U -> jac  = dmatrix (0, nsm1, 0, nrm1);
  U -> mass = dmatrix (0, nsm1, 0, nrm1);

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

/* ------------------------------------------------------------------------- *
 * make_*() --- Convert external to internal representation of a curve       *
 * ------------------------------------------------------------------------- */

struct line make_line (Element *elmt, Edge *edge)
{
  const int start = edge->start;
  const int skip  = edge->skip;
  const int n     = edge->np - 1;

  struct line c;

  c.pos[0].x = (*elmt->xmesh)[start];
  c.pos[0].y = (*elmt->ymesh)[start];

  c.pos[1].x = (*elmt->xmesh)[start + n*skip];
  c.pos[1].y = (*elmt->ymesh)[start + n*skip];
 
  return c;
}

struct arc make_arc (Element *elmt, Edge *edge, double radius)
{
  struct arc c;
  double *xi, alpha, theta, arclen;
  struct { double x, y; } p1, p2;
  
  p1.x   = (*elmt->xmesh)[edge->start];
  p1.y   = (*elmt->ymesh)[edge->start];
  p2.x   = (*elmt->xmesh)[edge->start + (edge->np-1) * edge->skip];
  p2.y   = (*elmt->ymesh)[edge->start + (edge->np-1) * edge->skip];
  arclen = sqrt(SQR(p2.x-p1.x) + SQR(p2.y-p1.y));

  if (fabs(alpha = radius/arclen) < 0.505)
    cubit_err ("arc radius is too small to fit");

  alpha = alpha > 0. ? -sqrt(SQR(alpha)-0.25) * arclen :
                        sqrt(SQR(alpha)-0.25) * arclen ;

  c.x      = (p1.x + p2.x) * .5 + alpha * (p2.y - p1.y) / arclen;
  c.y      = (p1.y + p2.y) * .5 + alpha * (p1.x - p2.x) / arclen;
  c.radius = fabs(radius);
  alpha    = fabs(alpha);

  p1.x    -= c.x;
  p1.y    -= c.y;

  c.theta.start = atan2 (p1.y, p1.x);
  c.theta.range = radius > 0. ? 2. * atan (.5 * arclen / alpha) :
                               -2. * atan (.5 * arclen / alpha) ;
  return c;
}

/* ------------------------------------------------------------------------- *
 * Bezier cubic splines.                                                     *
 *                                                                           *
 * spline_alloc (elmt, edge, pos)  -- Create a new spline given the coords   *
 *        of the intermediate control points, pos[] = { x1 y1 x2 y2 }.       *
 *                                                                           *
 * Bezier splines are specified by the values of a function at four nodes.   *
 *                                                                           *
 *                                                                           *
 *        f0          f1          f2           f3                            *
 *        o --------- o --------- o ---------- o                             *
 *        s=0         s=1/3       s=2/3        s=1                           *
 *                                                                           *
 * The spline itself is, of course, a cubic function that has four coeffi-   *
 * cients, f(s) = a s**3 + b s**2 + c s + d.  These coefficients are deter-  *
 * mines by requiring the spline to satisfy the following:                   *
 *                                                                           *
 *       f(0) = f0                                                           *
 *       f(1) = f3                                                           *
 *                                                                           *
 *      f'(0) = 3 (f1 - f0)                                                  *
 *      f'(1) = 3 (f3 - f2)                                                  *
 *                                                                           *
 * From these conditions, the value of the coefficients {a,b,c,d} can be     *
 * computed easily:                                                          *
 *                                                                           *
 *          d = f0                                                           *
 *          c = 3 (f1 - f0)                                                  *
 *          b = 3 (f2 - f1) - c                                              *
 *          a = f3 - f0 - c - b                                              *
 *                                                                           *
 * ------------------------------------------------------------------------- */

struct spline_t spline_alloc (Element *elmt, Edge *edge, double pos[])
{
  spline_t s;
  Point p0, p1, p2, p3;

  p0.x = (*elmt->xmesh)[edge->start];
  p0.y = (*elmt->ymesh)[edge->start];

  p1.x = pos[0];
  p1.y = pos[1];

  p2.x = pos[2];
  p2.y = pos[3];

  p3.x = (*elmt->xmesh)[edge->start + (edge->np-1)*edge->skip];
  p3.y = (*elmt->ymesh)[edge->start + (edge->np-1)*edge->skip];

  spline_compute(&s, p0, p1, p2, p3);

  { spline_t s1, s2;
    spline_cut (&s, &s1, &s2);
  }

  return s;
}

static void spline_compute 
(spline_t *s, Point p0, Point p1, Point p2, Point p3)
{
  s->x.d = p0.x;
  s->x.c = 3. * (p1.x - p0.x);
  s->x.b = 3. * (p2.x - p1.x) - s->x.c;
  s->x.a = p3.x - p0.x - s->x.c - s->x.b;

  s->y.d = p0.y;
  s->y.c = 3. * (p1.y - p0.y);
  s->y.b = 3. * (p2.y - p1.y) - s->y.c;
  s->y.a = p3.y - p0.y - s->y.c - s->y.b;
}
  
static void spline_getValue (spline_t *s, double t, double *x, double *y)
{
  *x = s->x.d + t*(s->x.c + t*(s->x.b + t*s->x.a));
  *y = s->y.d + t*(s->y.c + t*(s->y.b + t*s->y.a));
}

static void spline_fit (spline_t *s, int np, double *x, double *y)
{
  int i;
  double *xi; getops (np, &xi, 0, 0, 0);
  
  for (i = 0; i < np; i++) {
    double t = (xi[i] + 1.)/2.;
    spline_getValue (s, t, x + i, y + i);
  }

  return;
}

/* Get the value of a node (control point), node = 0 to 3 */

void spline_node (const spline_t *s, int node, double *x, double *y)
{
    double x0 = s->x.d;
    double x1 = x0  + s->x.c/3.;
    double x2 = x1 + (s->x.b + s->x.c)/3.;
    double x3 = s->x.a + s->x.b + s->x.c + s->x.d;

  switch (node) {
  case 0:
    *x = s->x.d;
    *y = s->y.d;
    break;
  case 1:
    *x = s->x.d + s->x.c/3.;
    *y = s->y.d + s->y.c/3.;
    break;
  case 2:
    *x = s->x.d + (2.*s->x.c + s->x.b)/3.;
    *y = s->y.d + (2.*s->y.c + s->y.b)/3.;
    break;
  case 3:
    *x = s->x.d + s->x.c + s->x.b + s->x.a;
    *y = s->y.d + s->y.c + s->y.b + s->y.a;
    break;
  default:
    cubit_err("spline_node: invalid node number");
    break;
  }
}

/* Cut a Bezier cubic spline into two equivalent splines over half the     *
 * interval.   The spline s1 is defined over (0,1/2), s2 is over (1/2,1).  */

void spline_cut (const spline_t *s, spline_t *s1, spline_t *s2)
{
  s1->x.a = s->x.a / 8.;
  s1->x.b = s->x.b / 4.;
  s1->x.c = s->x.c / 2.;
  s1->x.d = s->x.d;

  s1->y.a = s->y.a / 8.;
  s1->y.b = s->y.b / 4.;
  s1->y.c = s->y.c / 2.;
  s1->y.d = s->y.d;

  s2->x.a = s1->x.a;
  s2->x.b = 3.*s1->x.a + s1->x.b;
  s2->x.c = 3.*s1->x.a + 2.*s1->x.b + s1->x.c;
  s2->x.d = s1->x.a + s1->x.b + s1->x.c + s1->x.d;

  s2->y.a = s1->y.a;
  s2->y.b = 3.*s1->y.a + s1->y.b;
  s2->y.c = 3.*s1->y.a + 2.*s1->y.b + s1->y.c;
  s2->y.d = s1->y.a + s1->y.b + s1->y.c + s1->y.d;
}

/* ------------------------------------------------------------------------- */

static void lineto (int np, double *x, double *y, struct line line)
{
  int i;
  struct { double x, y; } slope, origin;
  double *xi;

  getops (np, &xi, 0, 0, 0);

  origin.x = line.pos[0].x;
  origin.y = line.pos[0].y;

  slope.x  = line.pos[1].x - line.pos[0].x;
  slope.y  = line.pos[1].y - line.pos[0].y;

  for (i = 0; i < np; i++) {
    x[i] = origin.x + slope.x * (xi[i] + 1.) * .5;
    y[i] = origin.y + slope.y * (xi[i] + 1.) * .5;
  }

  return;
}

static void arcto (int np, double *x, double *y, struct arc arc)
{
  int i;
  double *xi, theta;

  getops (np, &xi, 0, 0, 0);

  for (i = 0; i < np; i++) {
    theta = arc.theta.start + arc.theta.range * (xi[i] + 1.) * .5;
    x[i]  = arc.x + arc.radius * cos(theta);
    y[i]  = arc.y + arc.radius * sin(theta);
  }

  return;
}

/* ------------------------------------------------------------------------- */

#define MAX_NC 256

struct geofile_t geofile_alloc (const char *name)
{
  struct geofile_t g;
  char buf[BUFSIZ];
  FILE *fp;
  int i;

  memset(&g, 0, sizeof(geofile_t));
  g.name = strdup(name);

  if (!(fp = fopen(name,"r")))
    cubit_err ("unable to open the geometry file");

  /* Read the coordinates for this geometry */

  g.x = (double*) calloc (MAX_NC, sizeof(double));
  g.y = (double*) calloc (MAX_NC, sizeof(double));

  while (fgets(buf, BUFSIZ, fp)) 
    if (*buf != '#') break;
  i = 0;
  while (i <= MAX_NC && sscanf(buf, "%lf%lf", g.x+i, g.y+i) == 2) {
    ++i; if(!fgets(buf, BUFSIZ, fp)) break;
  }
  g.npts = i;

  /* Allocate space for the other things */

  g.sx     = (double*) calloc (g.npts, sizeof(double));
  g.sy     = (double*) calloc (g.npts, sizeof(double));
  g.arclen = (double*) calloc (g.npts, sizeof(double));

  /* Compute the actual arclength */

  g.arclen[0] = 0.;
  for (i = 0; i < g.npts-1; i++) {
    g.arclen[i+1] = g.arclen[i] +
      sqrt(SQR(g.x[i+1]-g.x[i]) + SQR(g.y[i+1]-g.y[i]));
  }

  /* Compute the spline information */

  spline (g.npts, 1.e30, 1.e30, g.arclen, g.x, g.sx);
  spline (g.npts, 1.e30, 1.e30, g.arclen, g.y, g.sy);

  fclose(fp);
  return g;
}
  
void _geofile_fit (geofile_t *g, Element *elmt, Edge *edge,
		   int np, double *xp, double *yp) {
  geofile_fit(g,elmt,edge,np,xp,yp);
}


static void geofile_fit (geofile_t *g, Element *elmt, Edge *edge,
			 int np, double *xp, double *yp)
{
  int i;

  Point p1 = Point_edge(elmt,  edge->id);
  Point a1 = Point_edge(elmt, (edge->id+3)%4);
  Point p2 = Point_edge(elmt, (edge->id+1)%4);
  Point a2 = Point_edge(elmt, (edge->id+2)%4);

  const int start = edge->start;
  const int skip  = edge->skip;

  double *eta = dvector(0,np);
  double *z; getops(np, &z, 0, 0, 0);

  /* Find the location where this edge intersects the curve */

  eta[0]    = geofile_search (g, a1, p1);
  eta[np-1] = geofile_search (g, a2, p2);

  /* Now generate the evaluation points */

  for (i = 1; i < np-1; i++)
    eta[i] = eta[0] + .5 * (eta[np-1] - eta[0]) * (z[i] + 1.);
  for (i = 0; i < np; i++) {
    xp[i] = splint(g->npts, eta[i], g->arclen, g->x, g->sx);
    yp[i] = splint(g->npts, eta[i], g->arclen, g->y, g->sy);
  }

  geofile_reset(g);
  free (eta);
  return;
}

/* ------------------------------------------------------------------------- */

static double geofile_search (geofile_t *g, Point a, Point p)
{
  const double tol = dparam("TOLCURV");
  const int    ip  = geofile_closest(g,p);

  double s[3], f[3];

  Vector ap;
  Vector_setValue (&ap, a, p);

  s[0] = g->arclen[ip];
  s[1] = g->arclen[ip+1];

  geofile_bracket (g, s, f, a, ap);
  if (fabs(f[1]) > tol)
    geofile_brent (g, s, a, ap, tol);

  return s[1];
}

static int geofile_closest (geofile_t *g, Point p)
{
  const double* x = g->x    + g->pos;
  const double* y = g->y    + g->pos;
  const int     n = g->npts - g->pos;

  double len[MAX_NC];
  int i;
  
  for (i = 0; i < n; i++)
    len[i] = sqrt(SQR(p.x-x[i]) + SQR(p.y-y[i]));

  i = idmin (n, len, 1) + g->pos;
  i = MIN (i, g->npts-2);

  /* If we found the same position and it's not the very first *
   * one, start the search over at the beginning again.  The   *
   * test for i > 0 makes sure we only do the recursion once.  */

  if (i && i == g->pos) { g->pos = 0; i = geofile_closest(g,p); }

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


static void geofile_bracket (geofile_t *g, double s[], double f[],
			     Point a, Vector ap)
{
  double ulim, u, r, q, fu;

  fa = geofile_angle (g, xa, a, ap);
  fb = geofile_angle (g, xb, a, ap);

  if (fb > fa) { SHFT (u, xa, xb, u); SHFT (fu, fb, fa, fu); }

  xc = xb + GOLD*(xb - xa);
  fc = geofile_angle (g, xc, a, ap);

  while (fb > fc) {
    r = (xb - xa) * (fb - fc);
    q = (xb - xc) * (fb - fa);
    u =  xb - ((xb - xc) * q - (xb - xa) * r) /
              (2.*SIGN(MAX(fabs(q-r),TINY),q-r));
    ulim = xb * GLIMIT * (xc - xb);

    if ((xb - u)*(u - xc) > 0.) {      /* Parabolic u is bewteen b and c */
      fu = geofile_angle(g, u, a, ap);
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
      fu = geofile_angle(g, u, a, ap);         /* default magnification. */

    } else if ((xc-u)*(u-ulim) > 0.) {   /* Parabolic fit is bewteen c   */
      fu = geofile_angle(g, u, a, ap);                     /* and ulim   */
      if (fu < fc) {
	SHFT  (xb, xc, u, xc + GOLD*(xc - xb));
	SHFT  (fb, fc, fu, geofile_angle(g, u, a, ap));
      }
    } else if ((u-ulim)*(ulim-xc) >= 0.) {  /* Limit parabolic u to the  */
      u   = ulim;                           /* maximum allowed value     */
      fu  = geofile_angle(g, u, a, ap);
    } else {                                       /* Reject parabolic u */
      u   = xc + GOLD * (xc - xb);
      fu  = geofile_angle(g, u, a, ap);
    }
    SHFT  (xa, xb, xc, u);      /* Eliminate the oldest point & continue */
    SHFT  (fa, fb, fc, fu);
  }
  return;
}

/* Brent's algorithm for parabolic minimization */

static double geofile_brent (geofile_t *g, double s[], Point ap, Vector app, 
			   double tol)
{
  int iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a  = MIN (xa, xc);               /* a and b must be in decending order */
  b  = MAX (xa, xc);
  d  = 1.;
  x  = w  = v  = xb;
  fw = fv = fx = geofile_angle (g, x, ap, app);

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
    fu = geofile_angle(g, u, ap, app);

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
  
  cubit_err ("too many iterations in brent()");
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

/* Compute the angle between the Vector ap and the vector connecting the     *
 * Point a to a particular position "s" along the curve.                     */

static double geofile_angle (geofile_t *g, double s, Point a, Vector ap)
{
  Point c = Point_new (splint(g->npts, s, g->arclen, g->x, g->sx),
		       splint(g->npts, s, g->arclen, g->y, g->sy));
  Vector ac;
  Vector_setValue (&ac, a, c);

  return 1. - ((ap.x * ac.x + ap.y * ac.y) / 
	       (Vector_length(ap) * Vector_length(ac)));
}

static int geofile_reset (geofile_t *g) { 
  return g->pos = 0; 
}
