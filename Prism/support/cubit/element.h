#ifndef ELEMENT_H
#define ELEMENT_H

/* ELEMENT
 * 
 * $Revision$
 *
 * Author:     R. D. Henderson
 *
 * This is the data structure used to represent a high-order (spectral)
 * element.  The geometry is curvilinear and described with an isoparametric
 * representation using the same high-order polynomial basis.  To allocate
 * an element you need to specify its order and the coordinates of its four 
 * corner nodes in clockwise order.  
 * ------------------------------------------------------------------------- */

#include "veclib/veclib.h"

#define ELEMENT_ID(elmt)     (elmt->id)
#define ELEMENT_KEY(elmt)    (elmt->key)
#define ELEMENT_FAMILY(elmt) (elmt->family)
#define ELEMENT_NR(elmt)     (elmt->nr)
#define ELEMENT_NS(elmt)     (elmt->ns)
#define ELEMENT_NZ(elmt)     (elmt->nz)
#define ELEMENT_NPTS(elmt)   (elmt->nr*elmt->ns*elmt->nz)

#define ELEMENT_MASS(elmt)   (elmt->mass[0])
#define ELEMENT_SIZE(elmt)   (dsum(elmt->nr*elmt->ns,*elmt->mass,1))
#define ELEMENT_DEPTH(elmt)  (Element_getLevel(elmt))
#define ELEMENT_REFINE(elmt) (elmt->workspace.ei[0])
#define ELEMENT_STATUS(elmt) (elmt->workspace.ei[1])
#define ELEMENT_DATA(elmt,u) (u->data[elmt->id])

#define ELEMENT_EDGE(elmt,i) (elmt->edge_list+i)
#define ELEMENT_VERT(elmt,i) (elmt->vert_list+i)

#define ELEMENT_ELIST(elmt)  (elmt->edge_list)
#define ELEMENT_VLIST(elmt)  (elmt->vert_list)
#define ELEMENT_NEXT(elmt)   (elmt->next)
#define ELEMENT_PREV(elmt)   (elmt->prev)

typedef struct element {            /* ....... ELEMENT Definition ........ */
  int              id           ;   /* id number                           */
  int              family       ;   /* family number from the base mesh    */
  int              key          ;   /* binary node key                     */
  int              nr, ns, nz   ;   /* resolution                          */

  int              *emap        ;   /* map from (nr,ns) to (nb,ni) storage */
  double           **xmesh      ;   /* mesh coordinates in physical space  */
  double           **ymesh      ;   /*                                     */

  double           **xr, **xs   ;   /* ----------------------------------- */
  double           **yr, **ys   ;   /*                                     */
  double           **jac        ;   /*       F A M I L Y   D A T A         */
  double           **rx, **ry   ;   /*             (Geometry)              */
  double           **sx, **sy   ;   /*                                     */
  double           **mass       ;   /* ----------------------------------- */

  union {                           /*     Element workspace (8 words)     */
    double ef [8];
    int    ei [8];
  } workspace;
  
  struct edge      *edge_list   ;   /* array of edges                      */
  struct vertex    *vert_list   ;   /* array of vertices                   */
  struct element   *next        ;   /* link to the next (if active)        */
  struct element   *prev        ;   /* link to the prev (if active)        */
} Element;

/* Prototypes */

Element* Element_alloc (int nr, int ns, int nz);
void     Element_free  (Element *elmt);

void Element_setGeometry (Element *elmt, double *xc, double *yc);
int  Element_setRefineOn (Element *elmt);
void Element_genxy       (Element *elmt);
void Element_dx          (Element *elmt, double *u, double *dx);
void Element_dy          (Element *elmt, double *u, double *dy);

/* Tools for interpolation */

void Element_project     (Element *elmt, const double *u, double *v, int quad);
void Element_project_2d  (Element *elmt, const double *u, double *v, int quad);

#endif





