#ifndef VDB_H
#define VDB_H

/* ------------------------------------------------------------------------- 
 * VDB: Voxel Database Class
 * 
 * $Revision$
 *
 * Author:      R. D. Henderson
 *
 * A VDB represents a collection of points in n-dimensional space (voxels).
 * Each voxel has a unique integer called its "key".  The VDB allows these
 * points to be treated as locations of shared memory, and provides a 
 * mechanism for combining local and remote data associated with the same
 * geometric position.
 *
 * References
 * ----------
 * R. D. Williams, "Voxel databases: a paradigm for parallelism with 
 * spatial structure," Concurrency 4(8), 1992.
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 * ------------------------------------------------------------------------- */

#include <stdio.h>

#define VDB_MAXDIM  3         /* maximum dimension */
#define VDB_SYNC    1000      /* message types */
#define VDB_COMBINE 1001      /* "    "    "   */

typedef struct vdb {
  int        dim;              /* spatial dimension */
  int        keyCount;         /* counter for available keys */
  int        keyMaximum;       /* maximum number of keys */
  int        keyHashTableSize; /* size of the key hash table */
  double     boxResolution;    /* resolution of the voxels */

  int        bufsiz;           /* communications buffer size */
  char*      buffer;           /* communications buffer */

  struct voxel **table;        /* hash table of voxel positions */
  struct alias **links;        /* table of remote voxel connections */
} vdb_t;

/* Create.  A vdb is created to hold a maximum number of voxels with a given *
 * resolution.  The parameter <boxsize> determines how close positions can   *
 * be before they are binned into the same voxel.                            */

vdb_t *vdb_alloc (int dim, int maxobjs, int bufsiz, double boxsize);
void   vdb_free  (vdb_t *vdb);

/* Position functions.  You can add a point (join), delete a point, or query *
 * the database to see if it already contains a position.  Join/delete/query *
 * all return the key for the position or -1 if it's not found.  After       *
 * making changes with join/delete the vdb must be updated by a call to      *
 * vdb_sync(), which returns the total number of stored positions.           */

int vdb_djoin   (vdb_t *vdb, const double coord[3]);
int vdb_ddelete (vdb_t *vdb, const double coord[3]);
int vdb_dquery  (const vdb_t *vdb, const double coord[3]);

int vdb_fjoin   (vdb_t *vdb, const float coord[3]);
int vdb_fdelete (vdb_t *vdb, const float coord[3]);
int vdb_fquery  (const vdb_t *vdb, const float coord[3]);

int vdb_ijoin   (vdb_t *vdb, const int coord[3]);
int vdb_idelete (vdb_t *vdb, const int coord[3]);
int vdb_iquery  (const vdb_t *vdb, const int coord[3]);

int vdb_sync    (vdb_t *vdb);

/* Communications.  The only form of communications handled by the vdb is   *
 * combining local/remote data.  Combining is the application of a commut-  *
 * ative operation (addition, multiplication, etc.) to a set of data.  You  *
 * must write a combine function to carry out the arithmetic and pass a     *
 * pointer to this function to vdb_combine().  Here is a simple example of  *
 * a CombineOp for adding two sets of doubles:                              *
 *                                                                          *
 *     void add (int nitems, size_t size, const void *p0, void *p1) {       *
 *       int n;                                                             *
 *       for (n = 0; n < nitems; n++)                                       *
 *         ((double*) p1)[n] += ((double*) p0)[n];                          *
 *     }                                                                    *
 *                                                                          *
 * Several common operations are pre-defined (see below).                   *
 *                                                                          *
 * The other communications operation is an "update" whereby data associat- *
 * ed with the secretary node for a group of points is distributed to all   *
 * other points in the group.                                               *
 *                                                                          *
 * Note that the array indicated by "ptr" should store "nitems" of the      *
 * given "size" at each array index, i.e. data for index "i" is stored at   *
 * the memory location (ptr + i*size*nitems).  The voxel keys are used to   *
 * manipulate chunks of data.                                               */

typedef void (*CombineOp)(size_t, int, const void*, void*);

int vdb_update  (const vdb_t *vdb, void *ptr, size_t size, int nitems);
int vdb_combine (const vdb_t *vdb, void *ptr, size_t size, int nitems, 
		 CombineOp op);

/* Some CombineOp's... [d = double, s = single, i = int], should be self-   *
 * explanatory. */

void comb_dsum  (size_t size, int nitems, const void *p0, void *p1);
void comb_dmult (size_t size, int nitems, const void *p0, void *p1);
void comb_isum  (size_t size, int nitems, const void *p0, void *p1);
void comb_imax  (size_t size, int nitems, const void *p0, void *p1);
void comb_imin  (size_t size, int nitems, const void *p0, void *p1);
void comb_iand  (size_t size, int nitems, const void *p0, void *p1);
void comb_ior   (size_t size, int nitems, const void *p0, void *p1);

/* Members & Rank.  vdb_members() fills an array with the hit count of each *
 * voxel.  This is the number of times vdb_add() was called for the same    *
 * position, minus the number of times vdb_delete() was called.  Return     *
 * value is the total number of positions in the database.                  *
 *                                                                          *
 * vdb_secr() fills an array with 1 or 0 depending on whether the local     *
 * voxel is designated as a "secretary" node.  Only one position in the     *
 * global set is designated as the secretary node.  This is used to imple-  *
 * ment a loop over voxels where each voxel must be visited only once, e.g. *
 * to compute a distributed dot product.                                    */

int vdb_members (const vdb_t *vdb, int *mult);
int vdb_secr    (const vdb_t *vdb, int *secr);

/* Info.  The function vdb_stats() prints a summary of the database cont-  *
 * ents (hash table efficiency, communications performance, etc.). The     *
 * function vdb_show() will dump the entire database onto an output stream *
 * for debugging or feedback.                                              */

int vdb_stats (const vdb_t *vdb, FILE *fp);
int vdb_show  (const vdb_t *vdb, FILE *fp);

/* Internal functions */

struct voxel *vdb_find (const vdb_t *vdb, const int coord[VDB_MAXDIM]);

/* ------------------------------------------------------------------------- *
 * The following provide compatibility with Roy William's VDB library        *
 * ------------------------------------------------------------------------- */

#define VDB struct vdb*

#define vdbcreate(dim,boxsize,_ig1,_ig2,maxobjs,_ig3,bufsize) \
     vdb_alloc(dim,maxobjs,bufsize,boxsize)
#define vdbdestroy(vdb) \
     vdb_free(vdb)

#define vdbdjoin(vdb,coord) \
     vdb_djoin(vdb,coord)
#define vdbfjoin(vdb,coord) \
     vdb_fjoin(vdb,coord)

#define vdbdpeek(vdb,coord) \
     vdb_dquery(vdb,coord)
#define vdbfpeek(vdb,coord) \
     vdb_fquery(vdb,coord)

#define vdbddelete(vdb,coord) \
     vdb_ddelete(vdb,coord)
#define vdbfdelete(vdb,coord) \
     vdb_fdelete(vdb,coord)

#define vdbsync(vdb) \
     vdb_sync(vdb)
#define vdbnmember(vdb,mult) \
     vdb_members(vdb,mult)
#define vdbsecr(vdb,secr) \
     vdb_secr(vdb,secr)

#define vdbcombine(vdb,ptr,op,size,nitems,stride) \
     vdb_combine(vdb,ptr,size,nitems,op)

#define vdbdsum comb_dsum
#define vdbisum comb_isum
#define vdbiand comb_iand
#define vdbimax comb_imax

#endif
