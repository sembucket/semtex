/*
 * Functions for dealing with patches
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
#include <math.h>
#include "mason.h"

/* Private functions */

static SegmentP 
  sort_segments (SegmentP list);
static PatchP
  sort_patches  (PatchP list);
static int 
  patch_cmp     (const void *p1, const void *p2);
static void
  segment_mesh  (SegmentP list, Patch *P),
  segment_link  (SegmentP masters, SegmentP slaves);

/* ------------------------------------------------------------------------ */

void linkPatches (Domain *omega)
{
  Patch *p;

  if (!omega->P) 
    return;
  else
    omega->P = sort_patches (omega->P);

  for (p = omega->P; p ; p = p->next) {

    if (!p->masters || !p->slaves) {
      sprintf   (error_buf, "Patch %d: missing master & slave edges", p->id);
      error_msg (error_buf);
    }

    p->masters = sort_segments (p->masters);     /* Sort by ID numbers */
    p->slaves  = sort_segments (p->slaves);
    p->origin  = p->masters->edge->right->vc;    /* Get the patch origin */

    segment_mesh (p->masters, p);                /* Compute meshes */
    segment_mesh (p->slaves , p);

    segment_link (p->masters, p->slaves);        /* Link the segments */
    segment_link (p->slaves, p->masters);
  }
  return;
}

Patch *makePatch (int id, Patch *link)
{
  Patch *new = (Patch*) calloc(1,sizeof(Patch));

  new->id   = id;
  new->next = link;

  return new;
}

/* Make a new segment */

Segment *makeSegment (int id, Segment *link, Edge *edge)
{
  Segment *new = (SegmentP) calloc(1,sizeof(Segment));

  new->id   = id;
  new->edge = edge;
  new->next = link;

  return new;
}

Patch *findPatch (int id, Domain *omega)
{
  Patch *p;

  for (p = omega->P; p ; p = p->next)
    if (p->id == id)
      break;

  return p;
}

Segment *findSegment (int id, Segment *slist)
{
  Segment *s;
  
  for (s = slist; s ; s = s->next)
    if (s->id == id)
      break;

  return s;
}

void showPatch   (Domain *omega)
{
  Patch *p = omega->P;

  while (p) {
    fprintf (stderr, "--------------------------------------\n"
	             "Patch ID %d, origin = (%g,%g)\n"
	             "--------------------------------------\n", 
	              p->id, p->origin.x, p->origin.y);

    fputs ("Masters:\n", stderr); showSegment (p->masters);
    fputs ("Slaves:\n" , stderr); showSegment (p->slaves);
    
    p = p->next;
  }
  return;
}
	    
void showSegment (Segment *s)
{
  int i, nb;

  while (s) {
    fprintf (stderr, "\tID = %d, (element,face) = (%d,%d)\n",
	     s->id, s->edge->iel, s->edge->id);
    fprintf (stderr, "\ts0 = %g, length = %g\n", s->s0, s->length);

    if ((nb = s->branches)) {
      fputs ("\tBranches to: ", stderr);
      for (i = 0; i < nb; i++)
	fprintf (stderr, "%d ", s-> branch_ID[i]);
      fputc ('\n', stderr);
    }

    if ((s = s->next))
      fputs ("\t.\n", stderr);
  }

  return;
}

/* ------------------------------------------------------------------------ *
 *                    P R I V A T E     F U N C T I O N S                   *
 * ------------------------------------------------------------------------ */

static PatchP sort_patches (Patch *list)
{
  int npatches = 0;
  int n;
  PatchP s, tmp;

  /* Count the patches in the list */

  for (s = list; s ; s = s->next) npatches++;            

  /* Create a new list for sorting */

  tmp    = (PatchP) malloc (npatches * sizeof(Patch));   
  for (n = 0, s = list; s ; s = s->next, n++) {
    memcpy (tmp+n, s, sizeof(Patch));
    free   (s);
  }

  qsort (tmp, npatches, sizeof(Patch), patch_cmp);

  for (n = 1; n < npatches; n++)
    tmp[n-1].next = tmp + n;

  tmp[--n].next = (PatchP) NULL;

  return tmp;
}

static SegmentP sort_segments (SegmentP list)
{
  int nsegs = 0;
  int n;
  SegmentP s, tmp;

  /* Count the segments in the list */

  for (s = list; s ; s = s->next) nsegs++;            

  /* Create a new list for sorting */

  tmp = (SegmentP) malloc (nsegs * sizeof(Segment));   
  s   = list;
  for (n = 0; n < nsegs; n++) {
    SegmentP next = s->next;
    s->next = NULL;
    memcpy (&tmp[n], s, sizeof(Segment));
    free (s);
    s = next;
  }

  qsort (tmp, nsegs, sizeof(Segment), patch_cmp);

  for (n = 0; n < nsegs-1; n++)
    tmp[n].next = &tmp[n+1];

  return tmp;
}

/* Compare the ID number of two patch-related structures */

static int patch_cmp (const void *p1, const void *p2)
{
  const int id1 = *((int*) p1);
  const int id2 = *((int*) p2);

  if      (id1 < id2)
    return -1;
  else if (id1 > id2)
    return  1;
  
  fprintf 
    (stderr, "warning: Patching ID numbers should be unique..."
             "check for %d's\n", id1);
  return 0;
}

/* Compute the mesh for each segment */

#define SQR(a) ((a)*(a))

static void segment_mesh (SegmentP list, Patch *P)
{
  double   x  = P->origin.x,
           y  = P->origin.y,
           s0 = 0.;
  SegmentP s  = list;
  Point    origin, endpoint;

  while (s) {
    if (s->edge->type == 'M') {
      origin   = s->edge->right->vc;
      endpoint = s->edge->left ->vc;
    } else {
      origin   = s->edge->left ->vc;
      endpoint = s->edge->right->vc;
    }      

    /* Compute the local segment offset and length */

    s->s0     = sqrt (SQR(origin.x-x) + SQR(origin.y-y)) + s0;
    s->length = sqrt (SQR(endpoint.x-origin.x) + SQR(endpoint.y-origin.y)); 

    /* Update these for the next segment */

    x   = endpoint.x;   
    y   = endpoint.y;
    s0 += s->length;
    s   = s->next;
  }  
  
  return;
}
    
#undef SQR

/* Link the "from" segments to the "to" segment */

#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )

#define MAX_BRANCHES 64
#define TOL_BRANCHES 1.e-6

static void segment_link (SegmentP from, SegmentP to)
{
  int      id_tmp[MAX_BRANCHES];
  double   lm, ls, len;
  SegmentP m, s;

  for (m = from; m ; m = m->next) {
    for (s = to; s  && m->branches < MAX_BRANCHES; s = s->next) {
      lm   = m->s0 + m->length;
      ls   = s->s0 + s->length;
      len  = MIN (lm, ls) - MAX (m->s0, s->s0);

      if (len < TOL_BRANCHES) continue;
	
      id_tmp [m->branches++] = s->id;
    }

    if (m->branches == MAX_BRANCHES)
      error_msg ("too many branches...please adjust and re-compile");

    if (m->branches)
      memcpy (m->branch_ID = (int*) malloc(m->branches * sizeof(int)),
	      id_tmp, m->branches * sizeof(int));
  }
  
  return;
}

#undef MIN
#undef MAX
#undef MAX_BRANCHES
#undef TOL_BRANCHES
