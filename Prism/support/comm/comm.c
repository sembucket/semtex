/*
 * COMM: Communications for Prism
 *
 * $Revision$
 *
 * Author:    R. D. Henderson
 *
 *
 * This file is just a stub for the communications library.  It includes the
 * header file "comm.h" and then grabs one of several drivers depending on
 * what's defined in CPPFLAGS.  To compile a "real" communications library
 * you need to define the token PARALLEL and one of the following:
 *
 *    MPI         Message Passing Interface
 *
 *    NX          Intel's NX message passing library
 *
 * Otherwise, the library compiles down to a set of dummy routines that do
 * not perform any communications.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "comm/comm.h"
#include "comm/config.h"

static int npes = 1;
static int rank = 0;

#if !defined(PARALLEL)
#  include "comm/driver.serial"
#else 
#
#  if defined(MPI)
#     include "comm/driver.mpi"
#  endif
#
#  if defined(NX)
#     include "comm/driver.nx"
#  endif 
#
#endif

/* Perform a ring test */

int comm_ringtest()
{
  int size = comm_size();
  int pid  = comm_rank();
  int i;

  const int type = 1000;

  int token = 1;
  int buf;

  comm_sync();

  if (pid == 0) {
    fprintf (stderr, "comm: starting ring test\n");
    comm_send (type, &token, sizeof(token), pid + 1);
    comm_recv (type, &buf, sizeof(token));
  } else {
    comm_recv (type, &buf, sizeof(token));
    fprintf (stderr, "comm: recived on node %d\n", pid);
    comm_send (type, &token, sizeof(token), (pid+1)%size);
  }

  comm_sync();

  if (pid == 0)
    fprintf (stderr, "comm: done\n");

  return 0;
}
