#ifndef ELEMENT_H
#define ELEMENT_H

/* Element
 *
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 *
 * Type declarations for an "element".
 * ----------------------------------------------------------------------- */

#include "speclib/edge.h"

#define ELEMENT_ID(elmt) (elmt->id)
#define ELEMENT_NR(elmt) (elmt->nr)
#define ELEMENT_NS(elmt) (elmt->ns)
#define ELEMENT_NZ(elmt) (elmt->nz)

#define ELEMENT_EDGE(elmt,i) (elmt->edges[i])
#define ELEMENT_XMESH(elmt)  (elmt->xmesh[0])
#define ELEMENT_YMESH(elmt)  (elmt->ymesh[0])

typedef struct element {            /* ....... ELEMENT Definition ........ */
  int              id           ;   /* ID (element) number                 */
  int              type         ;   /* Field type                          */
  int              nr, ns, nz   ;   /* Number of points in each direction  */
  int              frame        ;   /* Current frame number                */
  int              *emap        ;   /* Map from (nr,ns) to (nb,ni) storage */
  double           **base       ;   /* Pointer to the data storage area    */
  double           **field      ;   /*    (base) and the current frame     */
  double           **xmesh      ;   /* Mesh coordinates in physical space  */
  double           **ymesh      ;   /*                                     */

  double           **xr, **xs   ;   /* ----------------------------------- */
  double           **yr, **ys   ;   /*                                     */
  double           **jac        ;   /*       F A M I L Y   D A T A         */
  double           **rx, **ry   ;   /*             (Geometry)              */
  double           **sx, **sy   ;   /*                                     */
  double           **mass       ;   /* ----------------------------------- */

  struct edge      *edges       ;   /* Array of edges                      */
  struct element   *next        ;   /* Pointer to the next element         */
} Element;

/* Protypes */

void Element_grad (const Element *u, Element *dx, Element *dy);
void Element_dx   (const Element *u, Element *dx);
void Element_dy   (const Element *u, Element *dy);
  
#endif
