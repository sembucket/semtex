#ifndef EDGE_H
#define EDGE_H

/* EDGE
 *
 * $Revision$
 *
 * Author:     R. D. Henderson
 *
 * This is the data structure used to represent the edge of a spectral 
 * element.  It contains stride information to move along the edge, unit 
 * normals, and the quadrature information needed to compute surface 
 * integrals.
 * ------------------------------------------------------------------------- */

#define EDGE_ID(edge)     (edge->id)
#define EDGE_KEY(edge)    (edge->key)
#define EDGE_NP(edge)     (edge->np)
#define EDGE_NX(edge)     (edge->unx)
#define EDGE_NY(edge)     (edge->uny)
#define EDGE_DS(edge)     (edge->area)
#define EDGE_START(edge)  (edge->start)
#define EDGE_SKIP(edge)   (edge->skip)

typedef struct edge {               /* ........ EDGE definition .......... */
  int             id         ;      /* ID (face) number = [0..3]           */
  int             key        ;      /* Key to the global edge database     */
  int             dir        ;      /* Directional sense of the edge       */
  int             np         ;      /* Number of points along this edge    */
  double          pos[2]     ;      /* Position of the midpoint            */

  int             offset     ;      /* Starting boundary index offset      */
  int             start, skip;      /* Stride parameters                   */
  double          *unx, *uny ;      /* Unit outward normals                */
  double          *area      ;      /* Area associated with edge nodes     */

  struct bc       *bc        ;      /* Pointer to a boundary condition     */
  struct curve    *curve     ;      /* Pointer to a curve definition       */
  struct alias    *alias     ;      /* Pointer to an aliased edge          */
  struct vertex   *a, *b     ;      /* Pointer to two vertices             */
  struct edge     *next      ;      /* Pointer to the next edge            */
} Edge;

/* Prototypes */

Edge *Edge_valloc (int np, struct vertex *vert);
void  Edge_vfree  (Edge *vedge);

#endif
