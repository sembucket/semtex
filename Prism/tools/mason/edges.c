/*
 * Functions for operating on Edges
 *
 * RCS Information
 * ------------------------
 * $Author$
 * $Date$
 * $RCSfile$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mason.h"

/* Private Functions */

static void verify (Edge *edge, Domain *omega);

/* ------------------------------------------------------------------------ *
 * linkEdges() -- Link edges to vertices                                    *
 *                                                                          *
 * This is actually an initialization routine for the edges.  It just links *
 * each edge to its right/left nodes.  The sense of "left" and "right" is   *
 * defined by looking at the edges from the center of the element.  The     *
 * right node has the same local id as the edge, while the left node has    *
 * local id = edge_id + 1 % U->edges.                                       *
 *                                                                          *
 * The nodes along the edge of an element are shared by neighbors.  The     *
 * lower-numbered elements do the allocation, and the higher-numbered ele-  *
 * ments share a pointer the same memory location.                          *
 * ------------------------------------------------------------------------ */

void linkEdges (Domain *omega)
{
  const int tri = omega->tri;

  Element *U;
  Edge    *elist;
  VertexP *vlist;
  int     f;
  int     no_edges;

  for (U = omega->U; U ; U = U->next) {

    no_edges = U->edges;
    elist    = U->elist;
    vlist    = U->vlist;
    
    for (f = 0; f < no_edges; f++) {

      elist[f].next  = elist + f + 1;
      elist[f].right = vlist[f];
      elist[f].left  = vlist[(f + 1) % no_edges];

      /* Do we need to allocate edge pointers? */

      if (strchr("EP", elist[f].type)) {
	
	const int to_iel  = elist[f].bc.c.element;
	const int to_face = elist[f].bc.c.face;
	const int np      = elist[f].np;

	if (to_iel > elist[f].iel || 
	    (to_iel == elist[f].iel && to_face > elist[f].id)) {
	  elist[f].nodes = (Node*) calloc (np, sizeof(Node));
	  elist[f].dir   =  1;
	} else if (tri) {
	  elist[f].nodes = findEdge(omega,to_iel,to_face)->nodes;
	  elist[f].dir   =  1;
        } else {
	  elist[f].nodes = findEdge(omega,to_iel,to_face)->nodes + np-1;
	  elist[f].dir   = -1;
	}
      }

      /* Allocate a normal node array for non-connected edges */

      else {
	elist[f].nodes = (Node*) calloc (elist[f].np, sizeof(Node));
	elist[f].dir   = 1;
      }
    }
    elist[--f].next  = (Edge *) NULL;    /* Disconnect the last edge */
  }

  /* Now do a second pass to verify the connections */

  for (U = omega->U; U ; U = U->next)
    for (f = 0, no_edges = U->edges; f < no_edges; f++)
      verify (U->elist + f, omega);
  
  return;
}

/* ------------------------------------------------------------------------ *
 * findEdge() -- Search for a given edge in a mesh                          *
 *                                                                          *
 * Search for a given element/face combination in a mesh.  Returns a point- *
 * er to a valid edge or NULL if one can't be found.                        *
 * ------------------------------------------------------------------------ */

Edge *findEdge (Domain *omega, int element, int face)
{
  Element *U    = omega->U;
  Edge    *edge = (Edge *) NULL;

  while (U && U->id != element)            /* Start element search... */
    U    = U->next;

  if (U && face > 0 && face <= U->edges)   /* Check face number... */
    edge = U->elist + face - 1;

  if (!edge) {
    sprintf (error_buf, "edge %d of element %d not found\n", element, face);
    error_msg(error_buf);
  }

  return edge;
}

/* ------------------------------------------------------------------------ *
 * showEdges() -- Show the edge information for a domain                    *
 *                                                                          *
 * This is a debugging routine to display the connections for the edges     *
 * ------------------------------------------------------------------------ */

void showEdges (Domain *omega)
{
  Element *U;
  Edge    *edge;

  int n;

  fputs ("Showing edge information:\n", stderr);
  
  for (U = omega->U; U ; U = U->next) {
    fprintf (stderr, "Element %d\n", U->id);
    for (edge = U->elist; edge ; edge = edge->next) {
      const int dir = edge->dir;
      const int np  = edge->np;
      int       *p  = edge->nodes;

      fprintf (stderr, "%d - %c: ", edge->id, edge->type);
      fprintf (stderr, "%3d . ", *edge->right->node);
      for (n = 0; n < np; n++, p += dir)
	fprintf (stderr, "%3d ", *p);
      fprintf (stderr, ". %3d", *edge->left->node);

      switch (edge->type) {
      case 'E':
	fprintf (stderr, " --> Element %d, Edge %d", 
		 edge->bc.c.element, edge->bc.c.face);
	break;
      case 'P':
	fprintf (stderr, " <-> Element %d, Edge %d", 
		 edge->bc.c.element, edge->bc.c.face);
	break;
      case 'M': case 'S':
	fprintf (stderr, " *** Patched as %d.%d%c", 
		 edge->bc.p.patch, edge->bc.p.segment, tolower(edge->type));
	break;
      default:
	break;
      }

      putc ('\n',stderr);
    }
    if (U->next) fputs ("--------------------------------\n", stderr);
  }
  
  return;
}


/* ------------------------------------------------------------------------ *
 *                    P R I V A T E     F U N C T I O N S                   *
 * ------------------------------------------------------------------------ */

/* Check that element-element and periodic connections match up */

#define EPS 1.e-6

static void verify (Edge *edge, Domain *omega)
{
  const int type = edge->type;

  int  element, face;
  char buf [128];
  Edge *check;
  double dx1, dx2, dy1, dy2;

  if (type != 'E' && type != 'P') 
    return;

  /* Get the edge to check against */

  element = edge->bc.c.element;        
  face    = edge->bc.c.face;       
  check   = findEdge (omega, element, face);

  /* ... start the error checking ... */

  sprintf (buf, "element %d, face %d ", element, face);
  
  if (!check)
    error_msg (strcat (buf, "does not exist"));

  if (edge->type == 'E' && check->type != 'E')
    error_msg (strcat (buf, "is not an element-element boundary"));
  
  if (edge->type == 'P' && check->type != 'P')
    error_msg (strcat (buf, "is not a periodic boundary"));

  if (check->bc.c.element != edge->iel)
    error_msg (strcat (buf, "connects to the wrong element"));
  
  if (check->bc.c.face != edge->id)
    error_msg (strcat (buf, "connects to the wrong face"));

  /* Last check is for the geometry conditions */

  dx1 = fabs (edge->right->vc.x - check->left ->vc.x);
  dy1 = fabs (edge->right->vc.y - check->left ->vc.y);
  dx2 = fabs (edge->left ->vc.x - check->right->vc.x);
  dy2 = fabs (edge->left ->vc.y - check->right->vc.y);
     
  switch (type) {
  case 'E':
#if 0
    /* Check the difference in vertex coordinates */

    if (sqrt(dx1*dx1 + dy1*dy1) > EPS ||
	sqrt(dx2*dx2 + dy2*dy2) > EPS )
#else
    /* Check to make sure each edge points to the same vertex */

    if (edge->right != check->left || edge->left != check->right )
#endif
      error_msg (strcat (buf, "has the wrong vertex coordinates"));
    break;
  case 'P':
    if ((dx1 > EPS && dy1 > EPS) || (dx2 > EPS && dy2 > EPS))
      error_msg (strcat (buf, "differs in more than 1 coordinate direction"));
    break;
  default:
    break;
  }

  return;
}

#undef EPS

