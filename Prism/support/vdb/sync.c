/* 
 * vdb_t::sync implementation
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

/* Message buffer processing */

static int process (vdb_t *vdb, char *buf, int nvox, int done, status_t what)
{
  int n, p;

  const size_t header = 3*sizeof(int);
  const size_t packet = sizeof(voxel_t*) + 3*sizeof(int);

  /* First time through this loop we're transmitting the local buffer.  Each *
   * subsequent time we're passing another processor's buffer around the     *
   * ring.  The flag "spin" is turned off only when every buffer indicates   *
   * the local database has set done = 1.                                    */

  int spin = !done;
  int proc = comm_rank();

  /* Note: always send to the next processor in the ring */

  const int nprocs = comm_size();
  const int dest   = (proc + 1) % nprocs;

#ifdef DEBUG
  if (comm_rank()==0) {
    fprintf (stderr, "--------------------------\n");
    fprintf (stderr, "** -- ** %d %d\n", done, nvox);
  }
#endif

  for (p = 0; p < nprocs-1; p++) {
    const size_t nbytes = header + nvox*packet;
    char*        bufptr = buf;    

    comm_send (VDB_SYNC, buf, nbytes, dest);
    comm_recv (VDB_SYNC, buf, vdb->bufsiz);

    UPACK(bufptr, &nvox, sizeof(int));
    UPACK(bufptr, &done, sizeof(int)); 
    UPACK(bufptr, &proc, sizeof(int));

#ifdef DEBUG
    if (comm_rank()==0) {
      fprintf (stderr, "%2d -- %2d %d %d\n", p, proc, done, nvox);
    }
#endif

    for (n = 0; n < nvox; n++) {
      voxel_t *vox;   /* local entry */
      voxel_t *addr;  /* remote entry */
      int coord[3];   /* coordinates of remote entry */

      UPACK(bufptr, &addr, sizeof(voxel_t*));
      UPACK(bufptr, coord, sizeof(vox->coord));

      switch (what) {

      case VOXEL_DELETED:

	/* Look up the remote entry in the local database.  If it's present *
	 * then delete the existing alias to the remote entry.              */

	if (vox = vdb_find(vdb, coord))
	  voxel_unlink(vox, addr, proc);
	break;

      case VOXEL_NEW:
      case VOXEL_UPDATED:

	/* Look up the remote entry in the local database.  If it's present *
	 * then create an alias to the remote entry.                        */

	if (vox = vdb_find(vdb, coord))
	  voxel_link(vox, addr, proc);
	break;
	
      default:
	break;
      }
    }

    spin |= (!done);    /* ...see the note above */
  }

  return spin;
}

static int sync (vdb_t *vdb, status_t what) 
{
  const int    proc   = comm_rank();
  const int    size   = vdb->keyHashTableSize;
  const size_t packet = sizeof(voxel_t*) + 3*sizeof(int);
  const size_t nmax   = vdb->bufsiz - packet;

  int bin  = 0;
  int nvox = 0;
  int done = 0;
  int spin = 0;

  char *buf    = vdb->buffer;
  char *bufptr = buf;

  XPACK(bufptr, &nvox, sizeof(int));
  XPACK(bufptr, &done, sizeof(int));
  XPACK(bufptr, &proc, sizeof(int));

  for (bin = 0; bin < size; bin++) {
    voxel_t *vox = vdb->table[bin];
    if (!vox)
      continue;
    else do {

      if (what != VOXEL_ALL && what != vox->status)
	continue;

      /* Pack one voxel into the buffer for sending */
      XPACK(bufptr, &vox, sizeof(voxel_t*));
      XPACK(bufptr, vox->coord, sizeof(vox->coord));
      nvox++;

      /* Check for overflow on the communications buffer */
      if (bufptr - buf >= nmax) {
	bufptr = buf;
	XPACK(bufptr,&nvox,sizeof(int));
	XPACK(bufptr,&done,sizeof(int));
	XPACK(bufptr,&proc,sizeof(int));
	process(vdb, buf, nvox, done, what);
	nvox = 0;
      }
    } while (vox = vox->next);
  }

  /* Done.  Keep sending "done" message around the ring (spinning) until  *
   * everyone else is finished updating their own entries.                */

  done = spin = 1; 
  while (spin) {
    bufptr = buf;
    XPACK(bufptr,&nvox,sizeof(int));
    XPACK(bufptr,&done,sizeof(int));
    XPACK(bufptr,&proc,sizeof(int));

    spin = process(vdb, buf, nvox, done, what);
    nvox = 0;
  }

  return 0;
}

/* Remove all deleted voxels from the database */

static void unlink (vdb_t *vdb) 
{
  int bin;

  for (bin = 0; bin < vdb->keyHashTableSize; bin++) {
    voxel_t *vox  = vdb->table[bin];
    voxel_t *prev = vox;

    while (vox) {
      if (vox->status == VOXEL_DELETED) {
	voxel_t *next = vox->next;
	if (prev != vdb->table[bin])
	  prev->next = next;
	else
	  vdb->table[bin] = prev = next;
	voxel_free(vox);
	vox = next;
      } else {
	vox = vox->next;
      }
    }
  }
}


/* Rebuild the link tables for processor traversal and update the status of *
 * all entries to VOXEL_CONNECTED.                                          */

static void relink (vdb_t *vdb) 
{
  int bin;

  /* Clear the current list of aliases */

  memset (vdb->links, 0, comm_size()*sizeof(alias_t*));

  /* Loop through the set of voxels and construct a linked list connecting *
   * all of the aliases that reside on the same processor.   Update the    *
   * status of every voxel to VOXEL_CONNECTED.                             */

  for (bin = 0; bin < vdb->keyHashTableSize; bin++) {
    voxel_t *vox = vdb->table[bin];
    while (vox) {
      alias_t *alias = vox->link;
      while (alias) {
	const int p = alias->proc;
	alias->nextSameProc = vdb->links[p];
	vdb->links[p] = alias;
	alias = alias->next;
      }
      vox->status = VOXEL_CONNECTED;
      vox = vox->next;
    }
  }
}

/* Synchronize the database. */

int vdb_sync (vdb_t *vdb) 
{
  /* Disconnect old entries that have been deleted */
  sync(vdb, VOXEL_DELETED);

  /* Remove the deleted entries from the database */
  unlink(vdb);

  /* Connect new entries */
  sync(vdb, VOXEL_NEW);

  /* Connect old entries that have new links to the new links */
  sync(vdb, VOXEL_UPDATED);

  /* Rebuild the link tables for processor traversal */
  relink(vdb);
  
  return vdb->keyCount;
}
