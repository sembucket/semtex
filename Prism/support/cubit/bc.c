/*
 * BC implementation
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>

#include "bc.h"
#include "edge.h"
#include "element.h"

BC *BC_alloc (char type, int nitems)
{
  BC *bc     = (BC*) calloc(1,sizeof(BC));
  bc->type   = type;
  bc->nitems = nitems;

  return bc;
}

/* WARNING: This potentially leaks memory.  Since a BC does not allocate its *
 * own memory for info[] items, it does not de-allocate them either.  This   *
 * is probably no worse than a few bytes per BC, and only matters if the BC  *
 * stores expr strings or "other" stuff.                                     */

void BC_free (BC *bc) {
  free (bc);
}

/* Attach to the edge of an element */

void BC_attach (BC *bc, Element *elmt, Edge *edge)
{
  if (edge->bc)
    fprintf (stderr, "bc: edge is already attached!\n");
  else
    edge->bc = bc;
}
