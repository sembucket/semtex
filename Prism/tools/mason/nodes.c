/*
 * Functions for numbering the global nodes
 *
 * RCS Information
 * ------------------------
 * $Author$
 * $Date$
 * $RCSfile$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <string.h>
#include "mason.h"

static Node cross_link (Edge *edge, Domain *omega, Node pos);
static Node cleanup    (Domain *omega);

/* ------------------------------------------------------------------------ *
 * linkNodes() -- Global node numbering                                     *
 *                                                                          *
 * This function does the global node numbering.  The vertices are init-    *
 * tially unnumbered, and as this function passes through the mesh it sets  *
 * the values of each vertex and edge.  Element-element edges are numbered  *
 * only once by the first element; the later element just copy the nodes.   *
 *                                                                          *
 * Optimization can be turned on using the -O switch, or disables by -O0.   *
 * ------------------------------------------------------------------------ */
  
void linkNodes (Domain *omega)
{
  Element *U;

  for (U = omega->U; U ; U = U->next)       /* Preliminary numbering */
    linkNodesE (omega, U);

  if (0 < oplevel && oplevel < 3)           /* Optimize bandwidth */
    RCMop  (omega);
  else
    greedy (omega);

  omega->dof   = omega->nodes;
  omega->nodes = cleanup(omega);            /* Number delayed nodes */

  return;
}

/* Link the nodes in an element into the global mortar space */

void linkNodesE (Domain *omega, Element *U)
{
  Node pos = omega->nodes+1;
  Edge *edge;
  int  *p;

  for (edge = U->elist; edge ; edge = edge->next) {
    const int dir = edge->dir;
    int       np  = edge->np;

    if (!*(p = edge->right->node))
          *p = pos++;    

    switch (edge->type) {
    case 'S':
      if (!*(p = edge->nodes))
	while (np--) {
	  *p  = -1;            /* Flag for numbering later */
	   p += dir;
	}
      pos = cross_link (edge, omega, pos);
      break;

    default:
      if (!*(p = edge->nodes))         /* Has this edge been numbered ? */
	while (np--) {
	  *p  = pos++;
	   p += dir;
	}
      break;
    }

    if (!*(p = edge->left->node))
          *p = pos++;
  }
 
  omega->nodes = --pos;

  return;
}

/* ------------------------------------------------------------------------ *
 *                    P R I V A T E     F U N C T I O N S                   *
 * ------------------------------------------------------------------------ */

/* Number the delayed edges (slaves) */

static Node cleanup (Domain *omega)
{
  Element *U;
  Edge    *edge;
  Node     pos = omega->nodes+1, *p;
  
  for (U = omega->U; U; U = U->next)
    for (edge = U->elist; edge; edge = edge->next)

      if (*(p = edge->nodes) < 0) {
	int dir = edge->dir;
	int np  = edge->np;

	while (np--) {
	  *p  = pos++;
	   p += dir;
	}
      }

  return (omega->nodes = --pos);
}

/* Link across a slaved edge */

static Node cross_link (Edge *edge, Domain *omega, Node pos)
{
  PatchP   patch;
  SegmentP slave, master;
  int n, np, id, dir, *p;

  if (edge->type != 'S') {
    sprintf   (error_buf, "cross_link: edge %d is not a slave", edge->id);
    error_msg (error_buf);
  }

  if (!(patch = findPatch   (id = edge->bc.p.patch,   omega))) {
    sprintf   (error_buf, "cross_link: couldn't find Patch %d", id);
    error_msg (error_buf);
  }
  if (!(slave = findSegment (id = edge->bc.p.segment, patch->slaves))) {
    sprintf   (error_buf, "cross_link couldn't find slave %d.%ds", 
	                   edge->bc.p.patch, id);
    error_msg (error_buf);
  }

  for (n = 0; n < (*slave).branches; n++) {

    if (!(master = findSegment (id = slave->branch_ID[n], patch->masters))) {
      sprintf   (error_buf, "cross_link: failed searching for segment %d.%dm", 
		 patch->id, id);
      error_msg (error_buf);
    }
      
    dir = master->edge->dir;
    np  = master->edge->np;

    if (!*(p = master->edge->right->node))
          *p = pos++;

    if (!*(p = master->edge->nodes))
      while (np--) {
	*p  = pos++;
	 p += dir;
      } 

    if (!*(p = master->edge->left->node))
          *p = pos++;
  }

  return pos;
}
