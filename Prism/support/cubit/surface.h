#ifndef SURFACE_H
#define SURFACE_H

/* ------------------------------------------------------------------------- *
 * Surface                                                                   *
 *                                                                           *
 * Surfaces are implemented as a named list of (elmt,edge) pairs.  The list  *
 * is created as a linked list of these pairs, and a few generic list        *
 * operations are defined below.                                             *
 *                                                                           *
 * A loop over the segments in a surface is performed as follows:            *
 *                                                                           *
 *                                                                           *
 * s = Surface_alloc(mesh,"W");    // Group "W" = Walls                      *
 *                                                                           *
 * Surface_begin(s);                                                         *
 * while (!Surface_end(s)) {                                                 *
 *    Element *elmt = Surface_elmt(s);                                       *
 *    Edge    *edge = Surface_edge(s);                                       *
 *    ... Do some work ...                                                   *
 *    Surface_next(s);                                                       *
 * }                                                                         *
 *                                                                           *
 * Surface_free(s);                                                          *
 *                                                                           *
 * $Id$                                                                      *
 * ------------------------------------------------------------------------- */

#include "edge.h"
#include "element.h"
#include "mesh.h"

#define SURFACE_TYPE(s) (s->type)

typedef struct {
  char   type;
  struct segment *head;
  struct segment *curr;
} Surface;

/* Create and delete */

Surface *Surface_alloc(Mesh *mesh, char *typelist);
void     Surface_free (Surface *s);

/* Looping commands */

int Surface_begin (Surface *s);
int Surface_next  (Surface *s);
int Surface_more  (const Surface *s);
int Surface_end   (const Surface *s);
int Surface_len   (const Surface *s);

/* Extract the (elmt,edge) pair for current segment */

Element *Surface_elmt (const Surface *s);
Edge    *Surface_edge (const Surface *s);

#endif
