/*
 * Functions for processing the vertices
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
#include <math.h>
#include "mason.h"

/* Private functions */

static int     matchPoint  (Point p1, Point p2);

/* ------------------------------------------------------------------------ *
 * findVertex() -- Match a Vertex to a Point                                *
 *                                                                          *
 * This function matches an (x,y)-point to an existing member of a linked-  *
 * list of vertices.  It returns a pointer to the match or NULL of no match *
 * is found.                                                                *
 * ------------------------------------------------------------------------ */

Vertex *findVertex (Vertex *list, Point p)
{
  Vertex *v;

  for (v = list; v ; v = v->next)
    if (matchPoint (v->vc, p)) {
      v->mult++;
      break;
    }

  return v;
}

Vertex *makeVertex (Vertex *link, Point p)
{
  Vertex *new;

  new       = (Vertex *) calloc (1, sizeof(Vertex));
  new->id   = link ? link->id + 1 : 1;
  new->vc   = p;
  new->mult = 1;
  new->node = (Node *) calloc (1, sizeof(Node));
  new->next = link;

  return new;
}

/* ------------------------------------------------------------------------ *
 * linkVertex() -- Link the periodic vertices                               *
 *                                                                          *
 * This function resets the vertices linked across periodic boundaries.     *
 * Lower-numbered elements always take precedence over higher-numbered ele- *
 * ments.  When a node is reset, changes are propogated by modifying the    *
 * address of the shared node for that vertex.                              *
 * ------------------------------------------------------------------------ */

void linkVertex (Domain *omega)
{
  Element *U;
  Edge    *edge, *link;
  int      element, face;

  for (U = omega->U; U ; U = U->next) {
    for (edge = U->elist; edge ; edge = edge->next) {
      if (edge->type == 'P') {

	element =  edge->bc.c.element;
	face    =  edge->bc.c.face;

	if (element > edge->iel)
	  continue;
	else if (element == edge->iel && face > edge->id)
	  continue;

	/* We have to modify the vertices for this edge */

	if (!(link = findEdge (omega, element, face)))
	  error_msg ("search failed in linkVertex()");

	/* Reset the node pointers on this edge to match those  *
	 * on the "master" edge.  This propogates to every ele- *
	 * ment that shares these vertices.                     */

	edge->right->node = link-> left->node;
	edge-> left->node = link->right->node;
      }
    }
  }

  return;
}

/* ------------------------------------------------------------------------ *
 * showVertex() -- Display the list of vertices                             *
 *                                                                          *
 * This is a debugging routine to display the list of vertices and their    *
 * multiplicities.                                                          *
 * ------------------------------------------------------------------------ */

void showVertex (Domain *omega)
{
  Vertex *v;
  
  puts("Current vertex list:");
  for (v = omega->master; v ; v = v->next) {
    printf ("ID    = %d\n",  v->id);
    printf ("node  = %d\n", *v->node);
    printf ("(x,y) = %#14.7f %#14.7f\n", v->vc.x, v->vc.y);
    printf ("mult  = %d\n",  v->mult);
    
    if (v->next) puts ("----------------------");
  }

  return;
}

/* ------------------------------------------------------------------------ *
 *                    P R I V A T E     F U N C T I O N S                   *
 * ------------------------------------------------------------------------ */

/* Match two points */

#define EPS 1.e-6

static int matchPoint (Point p1, Point p2)
{
  int    ok = 0;
  double dx = p1.x - p2.x;
  double dy = p1.y - p2.y;

  if (sqrt (dx*dx + dy*dy) < EPS)
    ok = 1;

  return ok;
}

#undef EPS
