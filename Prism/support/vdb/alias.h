#ifndef ALIAS_H
#define ALIAS_H

/* ------------------------------------------------------------------------- 
 * alias_t
 *
 * $Id$
 *
 * Author: R. D. Henderson
 *
 * An alias is just a connection between a local voxel and a remote voxel. 
 * The data stored is the processor number and remote address of the link.
 * This way the keys on a remote processor can be renumbered but the con-
 * nections will remain valid.
 *
 * Copyright (c) 1998-1999 R. D. Henderson and Caltech.
 * ------------------------------------------------------------------------- */

typedef struct alias {         /* ---- voxel alias (remote connection) ---- */
  int           proc;          /* processor where the remote voxel lives    */
  struct voxel *remote;        /* address of the remote voxel               */
  struct voxel *local;         /* address of the local voxel                */
  struct alias *next;          /* pointer to the next in a linked list      */
  struct alias *nextSameProc;  /*   "    "    " w/same processor id         */
} alias_t;

/* Prototypes */

alias_t *alias_alloc();
void     alias_free (alias_t *alias);

#endif

