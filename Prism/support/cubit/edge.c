/* Edge implementation
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <assert.h>
#include "vertex.h"
#include "edge.h"

/* Edges can only be allocated in sets of 4 */

Edge *Edge_valloc (int np, Vertex *vert)
{
  Edge *vedge = (Edge*) calloc(4,sizeof(Edge));
  int id;

  assert(vedge);

  for (id = 0; id < 4; id++) {

    vedge[id].id   = id;
    vedge[id].np   = np;
    vedge[id].unx  = (double*) calloc(np,sizeof(double));
    vedge[id].uny  = (double*) calloc(np,sizeof(double));
    vedge[id].area = (double*) calloc(np,sizeof(double));
    vedge[id].a    = &vert[id];
    vedge[id].b    = &vert[(id+1)%4];

    switch (id) {
    case 0:
      vedge[id].offset = 0;
      vedge[id].start  = 0;
      vedge[id].skip   = 1;
      break;
    case 1:
      vedge[id].offset = np - 1;
      vedge[id].start  = np - 1;
      vedge[id].skip   = np;
      break;
    case 2:
      vedge[id].offset = np + np - 2;
      vedge[id].start  = np * np - 1;
      vedge[id].skip   = -1;
      break;
    case 3:
      vedge[id].offset = np * 2 + np - 3;
      vedge[id].start  = np * (np - 1);
      vedge[id].skip   = -np;
      break;
    }
  }

  for (id = 0; id < 3; id++)
    vedge[id].next = &vedge[id+1];

  return vedge;
}

void Edge_vfree (Edge *vedge)
{
  int id;

  for (id = 0; id < 4; id++) {
    free (vedge[id].unx);
    free (vedge[id].uny);
    free (vedge[id].area);
  }

  free (vedge);
}
