/*
 * Examples of the user-defined functions
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "cubit/cubit.h"
#include "mscope.h"

extern Domain Geometry;

/* Update Fields */  

int mscope_refine (Element *parent, Element *child[4])
{
  Mesh *mesh = Geometry.mesh;
  int i, k;

  /* Refine all solution fields */

  for (i = 0; i < Geometry.nfields; i++) {
    Field *u = Field_realloc(Geometry.solVector[i]);
    for (k = 0; k < 4; k++) {
      double *pdata = ELEMENT_DATA(parent,   u);
      double *cdata = ELEMENT_DATA(child[k], u);
      Element_project (parent, pdata, cdata, k);
    }
    Geometry.solVector[i] = u;
  }

  return 0;
}

/* ------------------------------------------------------------------------- *
 * This is an example of a user-defined BC handler.                          *
 *                                                                           *
 *                                                                           *
 * ------------------------------------------------------------------------- */

static void mscope_bc_vertex (BC *bc, Element *elmt, int id, double *data)
{
  const int  ip = id * (elmt->nr - 1);
  const int  nb = elmt->nr * 4 - 4;
  double xb[_MAX_NB], yb[_MAX_NB];
  
  dgathr (nb, *elmt->xmesh, elmt->emap, xb);
  dgathr (nb, *elmt->ymesh, elmt->emap, yb);

  switch (bc->type) {
  case 'V': case 'T': case 'F':
    *data = bc->info[0].value;
    break;
  case 'v': case 't': case 'f':
    scalar_set ("x", xb[ip]);
    scalar_set ("y", yb[ip]);
    *data = scalar (bc->info[0].expr);
    break;
  default:
    *data = 0.;
    break;
  }
}

static void mscope_bc_edge (BC *bc, Element *elmt, int id, double *data)
{
  const int ip = id * (elmt->nr - 1) + 1;
  const int nb = elmt->nr * 4 - 4;
  const int np = elmt->nr - 2;

  double xb[_MAX_NB];
  double yb[_MAX_NB];
  int i;

  dgathr (nb, *elmt->xmesh, elmt->emap, xb);
  dgathr (nb, *elmt->ymesh, elmt->emap, yb);

  switch (bc->type) {
  case 'V': case 'T': case 'F':
    for (i = 0; i < np; i++)
      data[i] = bc->info[0].value;
    break;
  case 'v': case 't': case 'f':
    vector_def ("x y", bc->info[0].expr);
    vector_set (np, xb + ip, yb + ip, data);
    break;
  default:
    for (i = 0; i < np; i++)
      data[i] = 0.;
    break;
  }
}

int mscope_bc (BC *bc, Element *elmt, int vertex, int edge, double *data)
{
  char type = bc->type;
  int  class;

  if (strchr("DVvTtW", type))
    class = DIRICHLET;
  else if (strchr("NFf", type))
    class = NEUMANN;
  else
    class = UNKNOWN;

  if (data != NULL) {
    if (-1 < vertex)
      mscope_bc_vertex (bc, elmt, vertex, data);
    if (-1 < edge)
      mscope_bc_edge   (bc, elmt, edge,   data);
  }
  
  return class;
}

