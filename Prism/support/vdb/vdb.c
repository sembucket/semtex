/*
 * vdb_t implementation
 *
 * $Revision$
 *
 * Author: R. D. Henderson
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "comm/comm.h"

#include "vdb.h"
#include "voxel.h"

/* ------------------------------------------------------------------------- */

/* Compute the smallest prime number greater than q */

static int prime (int q) 
{
  int p = (q < 20) ? 20 : q;
  while (--p > 19) {
    int sp = (int) sqrt((double) p), isp;
    for (isp = 2; isp<=sp; isp++)
      if ((p%isp) == 0)   /* p is not prime --- try (p-1) */
	break;            
    if (isp > sp) /* p is prime, so break out of the loop */
      break;
  }
  return p;
}

/* Hash function for voxel positions */

static int hash (const vdb_t *vdb, const int coord[VDB_MAXDIM]) {
#if VDB_MAXDIM==3
  int p = coord[0] + coord[1] + coord[2];
#else
  int p = 0;
  int i = 0;
  while (i < VDB_MAXDIM) p += coord[i++];
#endif
  return (p&0x7fffffff) % vdb->keyHashTableSize;
}

/* Lookup function */

static voxel_t *find (const vdb_t *vdb, const int coord[VDB_MAXDIM]) 
{
  const int bin = hash(vdb,coord);
  voxel_t  *vox = vdb->table[bin];

  while (vox) {
    if (memcmp(vox->coord, coord, sizeof(vox->coord)) == 0)
      break;
    vox = vox->next;
  }

  return vox;
}

/* Coordinate encoder(s) */

static void encode_d 
(const vdb_t *vdb, const double *xyz, int coord[VDB_MAXDIM])
{
  const int dim = vdb->dim;
  const double boxResolution = vdb->boxResolution;
  int i;

  for (i = 0; i < dim; i++)
    coord[i] = (int)(boxResolution*xyz[i]+0.5);
  while (i < VDB_MAXDIM)
    coord[i++] = 0;
}

static void encode_f
(const vdb_t *vdb, const float *xyz, int coord[VDB_MAXDIM])
{
  const int dim = vdb->dim;
  const float boxResolution = (float) vdb->boxResolution;
  int i;

  for (i = 0; i < dim; i++)
    coord[i] = (int)(boxResolution*xyz[i]+0.5);
  while (i < VDB_MAXDIM)
    coord[i++] = 0;
}

static void encode_i 
(const vdb_t *vdb, const int *xyz, int coord[VDB_MAXDIM]) 
{
  const int dim = vdb->dim;
  int i;

  for (i = 0; i < dim; i++)
    coord[i] = xyz[i];
  while (i < VDB_MAXDIM)
    coord[i++] = 0;
}

static int join (vdb_t *vdb, const int coord[VDB_MAXDIM]) 
{
  voxel_t *vox = find(vdb, coord);

  if (vox) {
    voxel_hit(vox);
  } else {
    const int bin = hash(vdb,coord);

    if (vdb->keyCount < vdb->keyMaximum) {
      vox       = voxel_alloc (coord);
      vox->key  = vdb->keyCount++;
      vox->next = vdb->table[bin];
      vdb->table[bin] = vox;
    } else {
      fprintf (stderr, "[%d]: no more keys: %d max\n", 
	       comm_rank(), vdb->keyMaximum);
      exit(-1);
    }
  }

  return vox->key;
}

static int query (const vdb_t *vdb, const int coord[VDB_MAXDIM]) 
{
  voxel_t *vox = find(vdb, coord);
  return vox ? vox->key : -1;
}

static int delete (vdb_t *vdb, const int coord[VDB_MAXDIM]) 
{
  voxel_t *vox = find(vdb, coord);

  if (vox) 
    return voxel_unhit(vox);
  else
    return -1;
}
  
/* ------------------------------------------------------------------------- */

vdb_t *vdb_alloc (int dim, int maxobjs, int bufsize, double boxsize)
{
  vdb_t *vdb = (vdb_t*) malloc(sizeof(vdb_t));

  /* default values for various parameters */
  if (bufsize==0) 
    bufsize = BUFSIZ;
  if (boxsize==0.0) 
    boxsize = FLT_EPSILON;
  
  /* data associated with the VDB */
  vdb->dim              = dim;
  vdb->keyCount         = 0;
  vdb->keyMaximum       = maxobjs;
  vdb->keyHashTableSize = prime(0.05*maxobjs);
  vdb->boxResolution    = 1./boxsize;

  /* communications buffer */
  vdb->bufsiz = bufsize;
  vdb->buffer = malloc(bufsize);

  /* hash table for the voxels and processor traversal table */
  vdb->table = (voxel_t**) calloc(vdb->keyHashTableSize, sizeof(voxel_t*));
  vdb->links = (alias_t**) calloc(comm_size(), sizeof(alias_t*));

  return vdb;
}

void vdb_free (vdb_t *vdb) 
{
  const int size = vdb->keyHashTableSize;

  int  bin;
  for (bin = 0; bin < size; bin++) {
    voxel_t *vox = vdb->table[bin];
    while (vox) {
      voxel_t *next = vox->next;
      voxel_free(vox);
      vox = next;
    }
  }

  free (vdb->buffer);
  free (vdb->table);
  free (vdb->links);
  free (vdb);
}

/* ------------------------------------------------------------------------- */

int vdb_djoin (vdb_t *vdb, const double *xyz) {
  int coord[VDB_MAXDIM];
  encode_d(vdb,xyz,coord);
  return  join(vdb,coord);
}

int vdb_dquery (const vdb_t *vdb, const double *xyz) {
  int coord[VDB_MAXDIM];
  encode_d(vdb,xyz,coord);
  return query(vdb,coord);
}

int vdb_ddelete (vdb_t *vdb, const double *xyz) {
  int coord[VDB_MAXDIM];
  encode_d (vdb,xyz,coord);
  return delete(vdb,coord);
}

/* ------------------------------------------------------------------------- */

int vdb_fjoin (vdb_t *vdb, const float *xyz) {
  int coord[VDB_MAXDIM];
  encode_f (vdb,xyz,coord);
  return   join(vdb,coord);
}

int vdb_fquery (const vdb_t *vdb, const float *xyz) {
  int coord[VDB_MAXDIM];
  encode_f(vdb,xyz,coord);
  return query(vdb,coord);
}

int vdb_fdelete (vdb_t *vdb, const float *xyz) {
  int coord[VDB_MAXDIM];
  encode_f (vdb,xyz,coord);
  return delete(vdb,coord);
}

/* ------------------------------------------------------------------------- */

int vdb_ijoin (vdb_t *vdb, const int *xyz) {
  int coord[VDB_MAXDIM];
  encode_i (vdb,xyz,coord);
  return   join(vdb,coord);
}

int vdb_iquery (const vdb_t *vdb, const int *xyz) {
  int coord[VDB_MAXDIM];
  encode_i(vdb,xyz,coord);
  return query(vdb,coord);
}

int vdb_idelete (vdb_t *vdb, const int *xyz) {
  int coord[VDB_MAXDIM];
  encode_i (vdb,xyz,coord);
  return delete(vdb,coord);
}

/* ------------------------------------------------------------------------- */

int vdb_members (const vdb_t *vdb, int *mult)
{
  if (mult) {
    const int nbins = vdb->keyHashTableSize;
    int bin;

    /* Store the local hit count in mult[] */

    for (bin = 0; bin < nbins; bin++) {
      voxel_t *vox = vdb->table[bin];
      while (vox) {
	mult[vox->key] = vox->hits;
	vox = vox->next;
      }
    }
    
    /* Combine across all processors */

    vdb_combine (vdb, mult, sizeof(int), 1, comb_isum);
  }

  return vdb->keyCount;
}

/* Select either the highest-ranking or lowest-ranking processor (semi- *
 * randomly) as the secretary node.                                     */ 

int vdb_secr (const vdb_t *vdb, int *secr)
{
  const int size = vdb->keyHashTableSize;

  int bin;
  int sum, m[2];

  voxel_t *vox;
  alias_t *alias;

  for (bin = 0; bin < size; bin++) {
    if (!(vox = vdb->table[bin]))
      continue;
    do
      secr[vox->key] = voxel_secr(vox);
    while 
	(vox = vox->next);
  }

  return 0;
}

/* Wrapper around the static "find" */

struct voxel *vdb_find (const vdb_t *vdb, const int coord[VDB_MAXDIM]) {
  return find(vdb,coord);
}
    
/* ------------------------------------------------------------------------- */

int vdb_stats (const vdb_t *vdb, FILE *fp)
{
  const int nprocs = comm_size();
  const int proc   = comm_rank();

  int size   = vdb->keyHashTableSize;
  int total  = 0;
  int empty  = 0;

  int nlinks_sum = 0, nlinks_max;
  int npacks_sum = 0, npacks_max;

  int ibuf[3];
  
  int  bin;
  for (bin = 0; bin < size; bin++) {
    voxel_t *vox = vdb->table[bin];
    if (vox) {
      while (vox) {
	total++;
	vox = vox->next;
      }
    } else
      empty++;
  }

  ibuf[0] = total; comm_isum (1, ibuf, &total);
  ibuf[0] = empty; comm_isum (1, ibuf, &empty);
  ibuf[0] = size;  comm_isum (1, ibuf, &size);

  for (bin = 0; bin < nprocs; bin++) {
    alias_t *alias = vdb->links[bin];
    if (alias) {
      nlinks_sum++;
      do
	npacks_sum++;
      while 
	(alias = alias->nextSameProc);
    }
  }

  ibuf[0] = nlinks_sum; comm_imax (1, ibuf, &nlinks_max);
  ibuf[0] = nlinks_sum; comm_isum (1, ibuf, &nlinks_sum);
  ibuf[0] = npacks_sum; comm_imax (1, ibuf, &npacks_max);
  ibuf[0] = npacks_sum; comm_isum (1, ibuf, &npacks_sum);

  if (proc == 0) {
    fprintf (fp, "vdb statistics:\n");
    fprintf (fp, "----------------------------------\n");
    fprintf (fp, "total entries   : %d\n", total);
    fprintf (fp, "hash table size : %d\n", size);
    fprintf (fp, "empty bins      : %4.2f%%\n", (100.*empty)/size);
    fprintf (fp, "avg occupancy   : %g\n", (double)total/(size-empty));
    fprintf (fp, "----------------------------------\n");
    fprintf (fp, "comm stats:  #links   #packets \n");
    fprintf (fp, "max          %6d   %6d\n", nlinks_max, npacks_max);
    fprintf (fp, "average      %6d   %6d\n", 
	     nlinks_sum/nprocs, npacks_sum/nprocs);
    fprintf (fp, "----------------------------------\n");
  }

  return 0;
}

int vdb_show (const vdb_t *vdb, FILE *fp)
{
  int bin;

  fprintf (fp, "vdb {\n");
  fprintf (fp, "\t dim %d\n", vdb->dim);
  fprintf (fp, "\t keyCount %d\n", vdb->keyCount);
  fprintf (fp, "\t keyMaximum %d\n", vdb->keyMaximum);
  fprintf (fp, "\t keyHashTableSize %d\n", vdb->keyHashTableSize);
  fprintf (fp, "\t boxResolution %g\n", vdb->boxResolution);
  fprintf (fp, "\t bufsize %d\n", vdb->bufsiz);
  
  for (bin = 0; bin < vdb->keyHashTableSize; bin++) {
    voxel_t *vox = vdb->table[bin];
    while (vox) {
      voxel_show (vox, fp);
      vox = vox->next;
    }
  }

  fprintf (fp, "}\n");
  return 0;
}
