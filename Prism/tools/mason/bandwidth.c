/* 
 * Bandwidth -- Compute the Mesh bandwidth
 *
 * This skips any node with an index < 0 (i.e. slaved edges and b.c.'s).
 * ---------------------------------------------------------------------- */

#include <stdio.h>
#include "mason.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* Private Functions */

static void minmax (int *min, int *max, int *pos, int np, int dir) 
{
  while (np--) {
    if (0 <= *pos) {
      *min = MIN(*min, *pos);
      *max = MAX(*max, *pos);
    }
    pos += dir;
  }
}

static void cross_bw (int *min, int *max, Domain *omega, Edge *edge)
{
  int b;

  PatchP   patch = findPatch   (edge->bc.p.patch,   omega);
  SegmentP slave = findSegment (edge->bc.p.segment, patch->slaves);

  for (b = 0; b < (*slave).branches; b++) {
    const int id     = slave->branch_ID[b];
    SegmentP  master = findSegment(id, patch->masters);

    if (!(master)) {
      sprintf   (error_buf, "cross_bw: failed searching for segment %d.%dm", 
		 patch->id, id);
      error_msg (error_buf);
    }
    
    minmax (min, max, master->edge->right->node, 1, 1);
    minmax (min, max, master->edge->left ->node, 1, 1);
    minmax (min, max, master->edge->nodes, 
	              master->edge->np, master->edge->dir);
  }
}

/* Compute the mesh bandwidth */

int bandwidth (Domain *omega)
{
  Element  *elmt;
  int i, bw = 0;
  
  for (elmt = omega->U; elmt; elmt = elmt->next) {
    int min = 0xffffff;
    int max = -min;

    for (i = 0; i < elmt->edges; i++) {
      minmax (&min, &max, elmt->elist[i].right->node, 1, 1);
      minmax (&min, &max, elmt->elist[i].left ->node, 1, 1);

      if (elmt->elist[i].type != 'S')
	minmax (&min, &max, 
		elmt->elist[i].nodes, elmt->elist[i].np, elmt->elist[i].dir);
      else
	cross_bw (&min, &max, omega, elmt->elist + i);
    }

    bw = MAX(bw, max-min);
  }

  return bw;
}

#undef MAX
#undef MIN
