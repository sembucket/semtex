/*
 * Vertex implementation
 *
 * ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <assert.h>
#include "vertex.h"

Vertex *Vertex_valloc()
{
  Vertex *vert = (Vertex*) calloc(4,sizeof(Vertex));
  int id;

  assert(vert);

  for (id = 0; id < 4; id++) 
    vert[id].id = id;
  for (id = 0; id < 3; id++)
    vert[id].next = &vert[id+1];

  return vert;
}

void Vertex_vfree (Vertex *vert) {
  free (vert);
}
