#ifndef VOXEL_H
#define VOXEL_H

/* ------------------------------------------------------------------------- 
 * voxel_t
 *
 * $Revision$
 *
 * Author:     R. D. Henderson
 *
 * Description:
 *
 * A voxel is a "volume element" of a 3-dimensional space.  Space is divided 
 * into discrete boxes (the voxels).  In this discrete space the voxels are
 * addressed using integer coordinates.
 *
 * Each voxel stores two integers: a "key" and a count of "hits".  The key is
 * a unique integer assigned to the voxel.  The count of hits tracks how many
 * times the position has been referencd.
 *
 * Also defined below is a structure called an "alias_t".  This is used to
 * set up a reference to a remote voxel, that is a voxel stored on another
 * processor with the same coordinates.
 *
 * Copyright (c) 1998-1999 R. D. Henderson and Caltech.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "alias.h"

typedef enum {              /* ----- voxel status types ----- */
  VOXEL_NEW        = 0,     /* ...newly created voxel         */
  VOXEL_UPDATED    = 1,     /* ...updated (alias added, etc.) */
  VOXEL_DELETED    = 2,     /* ...tagged for deletion         */
  VOXEL_CONNECTED  = 3,     /* ...connected                   */
  VOXEL_ALL        = 9      /* ...all status types            */
} status_t;

typedef struct voxel {      /* ----- voxel (volume element) ----- */
  int key;                  /* unique key for this voxel */
  int hits;                 /* number of hits */
  int coord[3];             /* coordinates of the voxel*/
  status_t status;          /* status of the voxel */
  struct alias *link;       /* pointer to the start of an alias list */
  struct voxel *next;       /* pointer to the next in a linked list */
} voxel_t;

/* Prototypes */

voxel_t *voxel_alloc (const int coord[3]);
void     voxel_free  (voxel_t *vox);

int voxel_hit    (voxel_t *vox);
int voxel_unhit  (voxel_t *vox);
int voxel_link   (voxel_t *vox, voxel_t *addr, int proc);
int voxel_unlink (voxel_t *vox, voxel_t *addr, int proc);
int voxel_secr   (const voxel_t *vox);
int voxel_show   (const voxel_t *vox, FILE *fp);

#endif
