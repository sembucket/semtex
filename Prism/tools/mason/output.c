/*
 * Functions for writing the output file
 *
 * RCS Information
 * ------------------------
 * $Author$
 * $Date$
 * $RCSfile$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include "mason.h"

/* ------------------------------------------------------------------------ *
 * output() -- Write the connectivity file                                  *
 *                                                                          *
 * This function writes the output file.  The format of this file is some   *
 * number of header lines containing useful information about the mesh      *
 * following by a listing of node numbers in the form:                      *
 *                                                                          *
 *           element  face   offset   node                                  *
 *                                                                          *
 * where "offset" is simply the distance along the edge.                    *
 * ------------------------------------------------------------------------ */

void output (Domain *omega, FILE *fp)
{
  Element *U;
  Edge    *edge;
  int i;

  fputs   (    "# Prism connectivity file\n"
	       "# \n", fp);
  fprintf (fp, "# Mesh bandwidth  = %d\n"
               "# Mesh DOF        = %d\n"  
               "# Mesh nodes      = %d\n",
	          bandwidth(omega), omega->dof, omega->nodes);

  if (oplevel)
  fprintf (fp, "# \n"
	       "# Optimization    = %d\n"
               "# First element   = %d\n", oplevel, opstart);

  fputs   (    "# \n"
	       "# elmt face offp node\n"
	       "# -------------------\n", fp);

  for (U = omega->U; U ; U = U->next)
    for (edge = U->elist; edge ; edge = edge->next) {
      const int dir = edge->dir;
      const int np  = edge->np;
      int       *p  = edge->nodes; 
      
      fprintf (fp, "  %4d %4d %4d %4d\n", 
	       U->id, edge->id, 1, *edge->right->node);
      for (i = 0; i < np; i++, p += dir)
	fprintf (fp, "  %4d %4d %4d %4d\n", 
		 U->id, edge->id, i+2, *p);
      fprintf (fp, "  %4d %4d %4d %4d\n", 
	       U->id, edge->id, i+2, *edge->left->node);
    }

  /* Patching information */

  if (omega->P) {
    Patch   *p;
    Segment *seg;

    fputs ("# Patching information \n"
	   "# \n"
	   "# type elmt face patch segment\n", fp);

    for (p = omega->P; p; p = p->next) {
      for (seg = p->masters; seg; seg = seg->next) {
	const Edge *edg = seg->edge;
	fprintf (fp, "     %c %4d %4d %4d %4d\n", 
		 edg->type, edg->iel, edg->id, 
		 edg->bc.p.patch, edg->bc.p.segment);
      }
      
      for (seg = p->slaves; seg; seg = seg->next) {
	const Edge *edg = seg->edge;
	fprintf (fp, "     %c %4d %4d %4d %4d\n", 
		 edg->type, edg->iel, edg->id, 
		 edg->bc.p.patch, edg->bc.p.segment);
      }
    }
  }

  fclose (fp);
}
