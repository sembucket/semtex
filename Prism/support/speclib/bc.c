/*
 * Boundary Conditions
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * -------------------------------------------------------------- */

#include <ctype.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "speclib/tree.h"

static Tree *BC_tags = NULL;

/* ------------------------------------------------------------------------- */

int BC_init(void)
{
  static char *Dirichlet_tags = "VvWDo";
  static char *Neumann_tags   = "AFfINO";
  char *p;

  if (BC_tags == NULL)
    BC_tags = create_tree(NULL, free);
  
  for (p = Dirichlet_tags; *p; p++)
    BC_defType(*p, DIRICHLET);
  for (p = Neumann_tags; *p; p++)
    BC_defType(*p, NEUMANN);

  return 0;
}


/* ------------------------------------------------------------------------- *
 * BC_free() -- free the memory allocated to boundary conditions             *
 *                                                                           * 
 * Although BCs are stored on what looks like a linked list, it's actually a * 
 * single contiguous array.  To free a list of BCs we first free up the      * 
 * memory allocated to storing boundary values, and then free the array of   * 
 * structures.                                                               * 
 * ------------------------------------------------------------------------- */

void BC_free (Bedge *list) 
{
  Bedge *bc = list;
  while (bc) {
    free(bc->bc.value);      /* extra memory allocated to this bc */
    bc = bc->next;
  }
  free(list);
}    

/* ------------------------------------------------------------------------- */

/* Count the number of boundary conditions */

int BC_count (Bedge *list) 
{
  int n = 0;
  while (list) { 
    n++; list = list->next; 
  }
  return n;
}

/* Duplicate a list of boundary conditions */

Bedge *BC_dup (Bedge *b1) 
{
  const int nbcs = BC_count(b1);

  Bedge *b2 = (Bedge*) calloc (nbcs, sizeof(Bedge));
  Bedge *bc = b1;
  int   ibc = 0;

  while (bc) {
    const size_t nbytes = (bc->edge->np)*(bc->elmt->nz)*sizeof(double);

    memcpy (b2 + ibc, bc, sizeof(Bedge));
    b2[ibc].bc.value = (double*) malloc (nbytes);
    memcpy (b2[ibc].bc.value, bc->bc.value, nbytes);

    ibc++; bc = bc->next;
  }

  /* Reset links */

  for (ibc = 0; ibc < (nbcs-1); ibc++)
    b2[ibc].next = &b2[ibc+1];

  return b2;
}

/* Find the boundaries of a given type in a list of bc's */

Bedge *BC_get (Bedge *list, char type)
{
  int    nbcs = 0;
  Bedge* sub  = NULL;
  Bedge* bc;
  
  /* Count the number of elements in the new list */

  for (bc = list; bc ; bc = bc->next)
    if (bc->type == type)
      nbcs++;

  /* Allocate memory for the new list */

  if (nbcs > 0) {
    sub  = (Bedge*) calloc (nbcs, sizeof(Bedge));
    bc   = list;
    nbcs = 0;
    
    /* Now copy them to the new list */

    do 
      if (bc->type == type) {
	memcpy (sub + nbcs, bc, sizeof(Bedge));
	sub [nbcs].next = sub + nbcs + 1;
	nbcs++;
      }
    while ((bc = bc->next));
    
    sub[nbcs-1].next = (Bedge*) NULL;
  }

  return sub;
}

/* 
 * Apply Dirichlet boundary conditions and check for singularity
 */

void BC_set (Bedge *Ubc, Field *U, BSystem *B)
{
  const int nb   = (U->nr + U->ns - 2) << 1;
  const int nel  = B->elements;
  const int bpts = B->bpts;
  int   singular = B->constant < FLT_EPSILON ? 1 : 0;
  int** bmap     = B->bmap;
  Bedge  *bc;

  int i, k;

  tempVector (u, bpts);

  /* This is a bit tricky.  To make sure boundary conditions are propogated
   * across the corner points, the boundary data is copied to the global
   * array "u" and then read back into every element, whether it contains
   * boundary edges or not.
   *
   * The alternative is to flag each element that touches a boundary.
   * ----------------------------------------------------------------------- */

  for (k = 0; k < nel; k++)                                 /* .. Initialize */
    for (i = 0; i < nb; i++)
      u [bmap[k][i]] = (*U[k].field) [U[k].emap[i]];

  for (bc = Ubc; bc; bc = bc->next) {                       /* ..... Correct */
    if (BC_getType(bc->type)==DIRICHLET) {
      if (isupper (bc->type) || bc->type == 'o') {
	int     np   = bc->edge->np;
	int    *node = bc->edge->bindex + bmap[bc->elmt->id];
	double *g    = bc->elmt->frame * np + bc->bc.value;

	singular = 0;  /* Turn off the singularity flag */

	for (i = 0; i < np; i++) 
	  u [ *(node++) ] = *g++;
      } 
    }
  }

  for (k = 0; k < nel; k++)                                 /* ..... Replace */
    for (i = 0; i < nb; i++)
      (*U[k].field) [U[k].emap[i]] = u [bmap[k][i]];

  B->singular = singular;

  freeVector (u);
}

/* "Learn" new boundary conditions from a field (all frames) */

void BC_learn (Bedge *Ubc, Field *U)
{
  const int nz   = U->nz;
  const int ntot = U->nr * U->ns * Field_count(U); 
  Bedge   *bedg;
  register int k;

  for (bedg = Ubc; bedg; bedg = bedg->next) {
    Edge    *edg  = bedg->edge;
    const int np  = edg ->np;
    const int id  = bedg->elmt->id;

    if (BC_getType(bedg->type)==DIRICHLET) {
      for (k = 0; k < nz; k++)
	ecopy (np, *U[id].base + ntot*k + edg->start, edg->skip,
	          bedg->bc.value + np*k, 1);
    }
  }

  return;
}

/* -------------------------------------------------------------------- *
 * BC_make() - Boundary edge structure setup                            *
 *                                                                      *
 * Boundary conditions are computed in physical space and then trans-   *
 * formed to Fourier space.                                             *
 *                                                                      *
 * This function no longer admits time-dependent boundary conditions.   *
 * -------------------------------------------------------------------- */

Bedge *BC_make (char bc, Element *elmt, Edge *edge, ...)
{
  va_list   ap;
  double    bv;
  char     *bf;
  Bedge    *new_bndry;

  /* Check to see if the boundary condition needs to be parsed */

  va_start (ap, edge);
  if (islower(bc)) {
    if ((bf = (strchr(va_arg(ap, char*), '='))) == (char*) NULL) {
      va_end(ap);
      speclib_error("illegal BC function definition:\n\t%s", bf);
    }
    while (isspace(*++bf));
  } else 
    bv = va_arg(ap, double);
  va_end(ap);

  new_bndry       = (Bedge*) malloc (sizeof(Bedge));
  new_bndry->type = toupper(bc);
  new_bndry->elmt = elmt;
  new_bndry->edge = edge;

  if (islower(bc)) {                                         /* expr */
    const int np  = edge->np;
    const int nz  = elmt->nz;
    double *value = dvector(0, np * nz - 1);
    double *x     = dvector(0, np);
    double *y     = dvector(0, np);
    int i;

    ecopy (np, *elmt->xmesh + edge->start, edge->skip, x, 1);
    ecopy (np, *elmt->ymesh + edge->start, edge->skip, y, 1);

    new_bndry->bc.value = value;
    
    if (nz == 1) {                                   /* 2-D version  */
      vector_def  ("x y", bf);
      vector_set  (np, x, y, value);
    } else {                                    
      int     NZ  = iparam("NZ");
      int     pid = option("procid");
      double *tmp = dvector(0, NZ);
      double *z   = zmesh  (NZ);

      /* 3-D Parallel    ||        Automatic Parallel Domain Mapping *
       *                                                             *
       * Just a reminder. NZ is the number of frames in the GLOBAL   *
       * domain, which must be a multiple of the # of processors.    *
       * The value U->nz is the number of LOCAL frames.  The b.c.'s  *
       * are computed, FFT'd, and then the local set of coefficients *
       * are stored.                                                 */

      for (i = 0; i < np; i++) {
	scalar_set  ("x", x[i]);
	scalar_set  ("y", y[i]);
	vector_def  ("z", bf);
	vector_set  (NZ, z, tmp);
	realft      (NZ/2, tmp, -1);
	dcopy       (nz, tmp + pid*nz, 1, value + i, np); 
      }
      free (tmp); 
      free (z);
    }
    free (x); 
    free (y);
  } else {                                           /* constant value */
    int    np  = edge->np;
    int    nz  = elmt->nz;
    int    pid = option ("procid");
    double *vv = (double*) calloc(np*nz, sizeof(double));

    /* BUG: This was previously defined to set the boundary condition
            if pid == 0 AND the patching flag was NOT set.             */

    if (pid == 0 || option("patch")) 
      dfill (np, bv, vv, 1);
    
    new_bndry->bc.value = vv;
  }

  return new_bndry;
}

/* ------------------------------------------------------------------------- */

BC_type BC_defType (char tag, BC_type type)
{
  char tagstr[2];
  Node *n;

  tagstr[0] = tag;
  tagstr[1] = '\0';

  if ((n = tree_search(BC_tags->root, tagstr))) {
    n->name  = strdup(tagstr);
    n->other = (void*) type;
  } else {
    n        = create_node(tagstr);
    n->other = (void*) type;
    tree_insert(BC_tags,n);
  }

  return type;
}

BC_type BC_getType (char tag)
{
  BC_type type = UNKNOWN;
  char tagstr[2];
  Node *n;

  tagstr[0] = tag;
  tagstr[1] = '\0';
  
  if ((n = tree_search(BC_tags->root, tagstr)))
    type = (BC_type) n->other;

  return type;
}


