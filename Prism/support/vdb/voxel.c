/*
 * voxel_t implementation
 *
 * $Revision$
 *
 * Author:       R. D. Henderson
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "voxel.h"
#include "alias.h"

static char *voxelStatus[] = {
  "VOXEL_NEW",
  "VOXEL_UPDATED",
  "VOXEL_DELETED",
  "VOXEL_CONNECTED"
};

/* ------------------------------------------------------------------------- */

voxel_t *voxel_alloc (const int coord[3]) 
{
  voxel_t *vox = (voxel_t*) malloc(sizeof(voxel_t));

  vox->key    = 0;
  vox->hits   = 1;
  vox->next   = NULL;
  vox->link   = NULL;
  vox->status = VOXEL_NEW;

  memcpy(vox->coord, coord, sizeof(vox->coord));
  return vox;
}

void voxel_free (voxel_t *vox) 
{
  alias_t *alias = vox->link;

  while (alias) {
    alias_t *next = alias->next;
    alias_free(alias);
    alias = next;
  }

  free (vox);
}

/* ------------------------------------------------------------------------- */

int voxel_link (voxel_t *vox, voxel_t *addr, int proc) 
{
  alias_t *alias;

  /* First check to see if we already have this one */

  for (alias = vox->link; alias; alias = alias->next) 
    if (alias->proc == proc && alias->remote == addr)
      return 0;

  /* Create a new alias */

  alias         = alias_alloc();
  alias->proc   = proc;
  alias->local  = vox;
  alias->remote = addr;
  alias->next   = vox->link;
  vox->link     = alias;

  /* Update the voxel status. A VOXEL_NEW doesn't need to change, but a *
   * VOXEL_CONNECTED must transition to VOXEL_UPDATED so that it can be *
   * shipped back out for linking with new remote voxels.               */

  if (vox->status == VOXEL_CONNECTED)
    vox->status = VOXEL_UPDATED;

  return (int) vox->status;
}

int voxel_unlink (voxel_t *vox, voxel_t *addr, int proc) 
{
  alias_t *alias = vox->link;
  alias_t *prev  = vox->link;

  while (alias) {

    if (alias->proc == proc && alias->remote == addr) {
      alias_t *next = alias->next;
      if (prev != vox->link)
	prev->next = next;
      else
	vox->link  = next;
      alias_free(alias);
      break;
    }

    prev  = alias;
    alias = alias->next;
  }

  return 0;
}

/* ------------------------------------------------------------------------- */

int voxel_compare (const voxel_t *v1, const voxel_t *v2) {
  return memcmp(v1->coord, v2->coord, sizeof(v2->coord));
}

int voxel_hit (voxel_t *vox) {
  return ++(vox->hits);
}

int voxel_unhit (voxel_t *vox) {
  if (vox->hits == 0)
    return 0;
  if (vox->hits == 1)
    vox->status = VOXEL_DELETED;
  return --(vox->hits);
}

/* Select either the highest-ranking or lowest-ranking processor (semi- *
 * randomly) as the secretary node.   Completely local computation.     */ 

int voxel_secr (const voxel_t *vox)
{
  const int proc = comm_rank();

  int secr = 1;    /* every voxel is a "secretary" by default */
  int sum, m[2];
  alias_t *alias;

  if (alias=vox->link) {
    m[0] = m[1] = sum = proc;
    do {
      if (alias->proc > m[0]) 
	m[0] = alias->proc;
      if (alias->proc < m[1]) 
	m[1] = alias->proc;
      sum += alias->proc;
    } while (alias=alias->next);

    if (m[sum%2] == proc)
      secr = 1;
    else
      secr = 0;
  }

  return secr;
}

/* ------------------------------------------------------------------------- */

int voxel_show (const voxel_t *vox, FILE *fp)
{
  fprintf (fp, "\t voxel {\n");
  fprintf (fp, "\t\t key %d\n", vox->key);
  fprintf (fp, "\t\t hits %d\n", vox->hits);
  fprintf (fp, "\t\t coord %d %d %d\n", 
	   vox->coord[0], vox->coord[1], vox->coord[2]);
  fprintf (fp, "\t\t status %s\n", voxelStatus[vox->status]);
  fprintf (fp, "\t\t addr %p\n", vox);

  if (vox->link) {
    alias_t *alias = vox->link;
    fprintf (fp, "\t\t links ");
    do
      fprintf (fp, "%d:%p ", alias->proc, alias->remote);
    while
      (alias = alias->next);
    fprintf (fp, "\n");
  }

  fprintf (fp, "\t }\n");
  return 0;
}


