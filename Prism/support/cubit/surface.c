/*
 * Surfaces
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>

#include "bc.h"
#include "edge.h"
#include "element.h"
#include "mesh.h"
#include "surface.h"

/* A segment of the surface */

struct segment {
  Element *elmt;
  Edge    *edge;
  struct segment *next;
}; 

/* ------------------------------------------------------------------------- */

Surface *Surface_alloc (Mesh *mesh, char *typelist)
{
  Surface *s = (Surface*) calloc(1,sizeof(Surface));

  Element *elmt;
  Edge    *edge;

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    for (edge = elmt->edge_list; edge; edge = edge->next) {
      if (edge->bc && strchr(typelist,edge->bc->type)) {
	struct segment *seg = (struct segment*) 
	  malloc(sizeof(struct segment));

	seg->elmt = elmt;
	seg->edge = edge;
	seg->next = s->head;
	s->head   = seg;
      }
    }
  }

  return s;
}

void Surface_free (Surface *s) 
{
  struct segment *seg = s->head;

  while (seg) {
    struct segment *next = seg->next;
    free(seg);
    seg = next;
  }

  free(s);
}

/* ------------------------------------------------------------------------- */

int Surface_begin (Surface *s) {
  s->curr = s->head;
  return 0;
}

int Surface_more (const Surface *s) {
  if (s->curr)
    return s->curr->next != NULL;
  else
    return 0;
}

int Surface_end (const Surface *s) {
  return s->curr == NULL;
}

int Surface_next (Surface *s) {
  int status = 0;

  if (s->curr)
    s->curr = s->curr->next;
  else
    status = -1;

  return status;
}

int Surface_len (const Surface *s) 
{
  int len = 0;
  struct segment *seg;

  for (seg = s->head; seg; seg = seg->next)
    len++;

  return len;
}

/* ------------------------------------------------------------------------- */

Element *Surface_elmt (const Surface *s) {
  return s->curr ? s->curr->elmt : NULL;
}

Edge *Surface_edge (const Surface *s) {
  return s->curr ? s->curr->edge : NULL;
}

    
