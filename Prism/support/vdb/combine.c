/*
 * vdb_t::combine
 *
 * $Revision$
 *
 * Author: R. D. Henderson
 * ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "comm/comm.h"

#include "vdb.h"
#include "voxel.h"

#define XPACK(buf,ptr,size) (memcpy(buf,ptr,size), buf += size)
#define UPACK(buf,ptr,size) (memcpy(ptr,buf,size), buf += size)

/* Ship out data to remote processors */

static void send (const vdb_t *vdb, char *ptr, size_t size, int nitems) 
{
  const size_t   data   = size*nitems;
  const size_t   packet = sizeof(voxel_t*) + data;
  const size_t   tail   = sizeof(voxel_t*) + sizeof(int);
  const size_t   nmax   = vdb->bufsiz - tail;

  const voxel_t *endtag = NULL;
  const int      proc   = comm_rank();
  const int      nprocs = comm_size();

  char *buf = vdb->buffer;

  int  dest;
  for (dest = 0; dest < nprocs; dest++) {
    alias_t *alias = vdb->links[dest];
    if (alias) {
      int done = 0;
      while (!done) {
	int   npack  = 0;
	char *bufptr = buf;
	
	do {
	  XPACK(bufptr,&(alias->remote),sizeof(voxel_t*));
	  XPACK(bufptr, ptr + alias->local->key*data, data);
	  npack++;
	} while 
	  ((alias=alias->nextSameProc) && ((bufptr-buf) < (nmax-packet)));
	  
	if (alias==NULL) done=1;

	XPACK(bufptr, &endtag, sizeof(endtag));
	XPACK(bufptr, &done,   sizeof(int));
	comm_send (VDB_COMBINE, buf, npack*packet+tail, dest);
      }
    }
  }
}

/* Receive remote data for combining */

static void recv 
(const vdb_t *vdb, char *ptr, size_t size, int nitems, CombineOp op)
{
  const size_t   data   = size*nitems;
  const size_t   tail   = sizeof(voxel_t*) + sizeof(int);

  const voxel_t *endtag = NULL;
  const int      proc   = comm_rank();
  const int      nprocs = comm_size();

  /* Note: There is an extra copy to unload data from the communications  *
   * buffer into the data-processing buffer "xpack".  The reason for this *
   * is to get around the potential alignment problem that occurs when    *
   * you pass an arbitrary address (buf+nbytes) into a CombineOp function *
   * that wants to read a particular data type from that address.  A copy *
   * to xpack works because malloc() returns an address that's valid for  *
   * ANY alignment type.  The hit for the extra copy should be negligible *
   * compared to the cost of getting the data here in the first place.    */

  char *buf   = vdb->buffer;
  char *xpack = (char*) malloc(data);

  int done   = 0;
  int nlinks = 0;
  int n;

  for (n = 0; n < nprocs; n++)
    if (vdb->links[n])
      nlinks++;

  while (done != nlinks) {
    int status = 0;
    char *bufptr = buf;
    voxel_t *addr;

    comm_recv (VDB_COMBINE, buf, vdb->bufsiz);

    for (UPACK(bufptr, &addr, sizeof(voxel_t*)); addr != endtag;
	 UPACK(bufptr, &addr, sizeof(voxel_t*)))
      {
	UPACK(bufptr, xpack, data);
	(*op)(size, nitems, xpack, ptr + addr->key*data);
      }

    UPACK(bufptr, &status, sizeof(int));
    done += status;
  }

  free (xpack);
}
  
int vdb_combine 
(const vdb_t *vdb, void *ptr, size_t size, int nitems, CombineOp op)
{
  send (vdb, (char*) ptr, size, nitems);
  recv (vdb, (char*) ptr, size, nitems, op);

  return comm_sync();
}

/* ------------------------------------------------------------------------- */

/* Add doubles */

void comb_dsum (size_t size, int nitems, const void *p0, void *p1) {
  int n;

  if (size != sizeof(double))
    return;

  for (n = 0; n < nitems; n++)
    ((double*) p1)[n] += ((double*) p0)[n];
}

/* Multiply doubles */

void comb_dmult (size_t size, int nitems, const void *p0, void *p1) {
  int n;

  if (size != sizeof(double))
    return;

  for (n = 0; n < nitems; n++)
    ((double*) p1)[n] *= ((double*) p0)[n];
}

/* Add integers */

void comb_isum (size_t size, int nitems, const void *p0, void *p1) {
  int n;

  if (size != sizeof(int))
    return;

  for (n = 0; n < nitems; n++)
    ((int*) p1)[n] += ((int*) p0)[n];
}

/* Logical AND */

void comb_iand (size_t size, int nitems, const void *p0, void *p1) {
  int n;

  if (size != sizeof(int))
    return;

  for (n = 0; n < nitems; n++)
    ((int*) p1)[n] &= ((int*) p0)[n];
}

/* Logical OR */

void comb_ior (size_t size, int nitems, const void *p0, void *p1) {
  int  n;

  if (size != sizeof(int))
    return;

  for (n = 0; n < nitems; n++)
    ((int*) p1)[n] |= ((int*) p0)[n];
}

/* Compute the integer minimum */

void comb_imin (size_t size, int nitems, const void *p0, void *p1) {
  int n;

  if (size != sizeof(int))
    return;

  for (n = 0; n < nitems; n++) {
    const int i0 = ((int*) p0)[n];
    const int i1 = ((int*) p1)[n];
    if (i0 < i1)
      ((int*) p1)[n] = i0;
  }
}

/* Compute the integer maximum */

void comb_imax (size_t size, int nitems, const void *p0, void *p1) {
  int n;

  if (size != sizeof(int))
    return;

  for (n = 0; n < nitems; n++) {
    const int i0 = ((int*) p0)[n];
    const int i1 = ((int*) p1)[n];
    if (i0 > i1)
      ((int*) p1)[n] = i0;
  }
}



