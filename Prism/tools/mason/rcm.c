/*
 * Bandwidth optimizer based on RCM
 * 
 * RCS Information
 * ------------------------
 * $Author$
 * $Date$
 * $RCSfile$
 * $Revision$
 * ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mason.h"

#if defined(CRAY)
#define genrcm_  GENRCM
#define rcmroot_ RCMROOT
#endif

#if defined(HPUX)
#define genrcm_  genrcm
#define rcmroot_ rcmroot
#endif

/* Internal variables */

typedef struct node {          /* ....... Node of the Graph ........ */
  int          n;
  struct node *next;
} node;

/* Build the node map for an element */

static int *build_map (Domain *omega, Element *elmt, int *npts)
{
  int n    = 0;
  int nmax = 64;
  int *map = (int*) malloc (nmax*sizeof(int));

  Edge  *edge;

  for (edge = elmt->elist; edge; edge = edge->next) {
    switch (edge->type) {
    case 'S': {
      Patch   *patch = findPatch  (edge->bc.p.patch, omega);
      Segment *slave = findSegment(edge->bc.p.segment, patch->slaves);
      const int nb   = slave->branches;
      register int i;

      for (i = 0; i < nb; i++) {
	Segment *master = findSegment (slave->branch_ID[i], patch->masters);

	int  dir = master->edge->dir;
	int  np  = master->edge->np;
	int *p   = master->edge->nodes;

	if (nmax - n < np + 2) {
	  nmax += 64;
	  map   = (int*) realloc (map, nmax*sizeof(int));
	}

	map[n++] = *master->edge->right->node;
	while (np--)
	  { map[n++] = *p; p += dir; }
	map[n++] = *master->edge->left->node;
      }
      break;
    }
    default: {
      int  dir = edge->dir;
      int  np  = edge->np;
      int *p   = edge->nodes;

      if (nmax - n < np + 2) {
	nmax += 64;
	map   = (int*) realloc (map, nmax*sizeof(int));
      }

      map[n++] = *edge->right->node;
      while (np--)
	{ map[n++] = *p; p += dir; }
      map[n++] = *edge->left->node;
    }
    }
  }

  *npts = n;
  return map;
}
 
/* Add node j to node i's adjacency list */

static void node_add (node **list, int i, int j)
{
  int found = 0;

  if (i != j) {
    node *head = list[i], *curr;

    for (curr = head; curr; curr = curr->next)
      if ((found = (j == (curr->n)))) break;
    
    if (!found) {
      curr       = (node*) calloc (1, sizeof(node));
      curr->n    = j;
      curr->next = head;
      list[i]    = curr;
    }
  }
}

static int node_count (node *list)
{
  int   n = 0;
  node *p = list;

  while (p) { n++; p = p->next; }
  
  return n;
}

static void adjncy_fill(node **list, int *adjncy, int *xadj, int size, int neq)
{
  int i, k;
  node  *p;

  k = 0; for (i = 0; i < neq; i++) {
    xadj[i] = k + 1;
    for (p = list[i]; p; p = p->next)
      adjncy[k++] = p->n + 1;
  }

  if (k != size) fprintf (stderr, "bad adjacency table\n");

  adjncy[k] = 0;
  xadj[neq] = k + 1;
}


static int adjncy_build (node **list, int neq, Domain *omega)
{
  Element *elmt;
  int      ip, np, size;
  register int i, j;

  for (elmt = omega->U; elmt; elmt = elmt->next) {
    int *map = build_map (omega, elmt, &np);

    for (i = 0; i < np; i++) {
      ip = map[i];
      for (j = 0; j < np; j++)
	node_add (list, ip, map[j]);
    }
    free (map);
  }

  size = 0;
  for (i = 0; i < neq; i++)
    size += node_count (list[i]);

#if DEBUG
  printf ("Adjacency list, size = %d:\n", size);
  for (i = 0; i < neq; i++) {
    node *p = list[i];
    printf ("%2d ", i);
    while (p) { printf ("%2d ", p->n); p = p->next; }
    putchar ('\n');
  }
#endif

  return size;
}

static void adjncy_free (node **list, int neq)
{
  node *n;
  int   i;

  for (i = 0; i < neq; i++) {
    while ((n = list[i])) { 
      list[i] = n->next; 
      free(n); 
    }
  }
}


static void map_load (Domain *omega, int *map)
{
  Element *elmt;
  Edge    *edge;
  Vertex  *vert;
  int     *p;

  for (vert = omega->master; vert; vert = vert->next)
    if (*(p = vert->node) > 0)
      vert->node = map + *p - 1;
  
  for (elmt = omega->U; elmt; elmt = elmt->next)
    for (edge = elmt->elist; edge; edge = edge->next) {
      if (*(p = edge->nodes) > 0)
	edge->nodes = map + *p - 1;
    }
}

static void map_invert (int neq, int *perm, int *map) {
  int  i;
  for (i = 0; i < neq; i++)
    map[perm[i]-1] = i+1;
}

/* ----------------------------------------------------------------------- */

void RCMop (Domain *omega)
{
  node **list;
  int  *adjncy, *xadj, *xls, *perm, *mask, *map, *root, size;
  int  i;   

  Element *elmt;
  int      bw, id;

  int  neq  = omega->nodes;
  int  npts = omega->nodes;

  /* External (Fortran) functions */

  void genrcm_();
  void rcmroot_();

  /* Default node numbering */

  map  = (int*) calloc (npts, sizeof(int));
  for (i = 0; i < npts; i++)
    map[i] = i;
  map_load (omega, map);

  /* Now all elements have pointers into "map" */

  list = (node**) calloc (npts, sizeof(node*));
  size = adjncy_build (list, npts, omega);

  /* Allocate memory for RCM routines */

  adjncy = (int*) calloc (size+1, sizeof(int));
  xadj   = (int*) calloc (neq+1, sizeof(int));
  perm   = (int*) calloc (neq, sizeof(int));
  mask   = (int*) calloc (neq, sizeof(int));
  xls    = (int*) calloc (neq, sizeof(int));
  root   = (int*) calloc (neq, sizeof(int));

  /* Number */

  adjncy_fill (list, adjncy, xadj, size, neq);
  genrcm_     (&neq, xadj, adjncy, perm, mask, xls);

  map_invert (neq, perm, map);
  bw = bandwidth(omega);
  id = -1;

#ifdef TEST
  printf ("Bandwidth [init] = %d\n", bw);
#endif

  /* Try each vertex of each boundary element as a starting node */

  if (oplevel > 1) {
    for (elmt = omega->U; elmt; elmt = elmt->next) {
      Edge *edge;
      for (edge = elmt->elist; edge; edge = edge->next)
	if (strchr("EPMS", edge->type) == NULL) {
	  root[*edge->right->node] = *edge->right->node + 1;
	  root[*edge->left ->node] = *edge->left ->node + 1;
	}
    }
    
    for (i = 0; i < neq; i++) {
      if (root[i]) {
	int bw_i;
	rcmroot_ (&neq, xadj, adjncy, perm, mask, xls, root+i);
	
	map_invert (neq, perm, map);
	if ((bw_i = bandwidth (omega)) < bw) {
	  bw = bw_i;
	  id = i;
	}
      }
    }

#ifdef TEST
    printf ("Bandwidth [best] = %d from node %d\n", bw, id > 0 ? root[id] : 0);
    exit (1);
#endif

    if (id >= 0) {
      opstart = root[id];
      rcmroot_ (&neq, xadj, adjncy, perm, mask, xls, &opstart);
      map_invert (neq, perm, map);
    }
  }

  adjncy_free (list, neq);

  free (adjncy);
  free (xadj);
  free (perm);
  free (mask);
  free (root);

  return;
}

