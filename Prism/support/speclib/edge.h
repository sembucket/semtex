#ifndef EDGE_H
#define EDGE_H

/* Edge
 *
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 *
 * ------------------------------------------------------------------------- */

#define EDGE_ID(edge)    (edge->id)
#define EDGE_NPTS(edge)  (edge->np)
#define EDGE_START(edge) (edge->start)
#define EDGE_SKIP(edge)  (edge->skip)
#define EDGE_NX(edge,i)  (edge->unx[(i)])
#define EDGE_NY(edge,i)  (edge->uny[(i)])
#define EDGE_DS(edge,i)  (edge->area[(i)])

typedef struct edge {               /* ........ EDGE definition .......... */
  int             id         ;      /* ID (face) number = [0..3]           */
  int             np         ;      /* Number of points along this edge    */
  int             bindex     ;      /* Starting boundary index => bmap     */
  int             start, skip;      /* For moving along edge of (*field)[] */
  double          *unx, *uny ;      /* Unit outward normals                */
  double          *area      ;      /* Area associated with edge nodes     */
  struct edge     *next      ;      /* Pointer to the next one             */
} Edge;

/* These are the methods to copy to/from an edge array and an element array */

void edge_gathr (const Edge *edge, const double *u, double *ue);
void edge_scatr (const Edge *edge, const double *ue, double *u);

/* These are various methods for performing line and surface integrals along *
 * an edge.  Each one takes an array of data (or two) and performs a quad-   *
 * rature to evaluate the integral.  The input data corresponds to:          *
 *                                                                           *
 *    _s     scalar function of position f(s)                                *
 *    _nx    scalar function of position u(s) dotted agains the unit normal  *
 *           in the x-direction                                              *
 *    _ny    scalar function of position u(s) dotted agains the unit normal  *
 *           in the y-direction                                              *
 *    _v     vector function of position i u(s) + j v(s); in this case the   *
 *           vector is dotted with the unit outward normal at each point.    */

double edge_integrate_s (const Edge *edge, const double *f);
double edge_integrate_nx(const Edge *edge, const double *u);
double edge_integrate_ny(const Edge *edge, const double *u);
double edge_integrate_v (const Edge *edge, const double *u, const double *v);

#endif
