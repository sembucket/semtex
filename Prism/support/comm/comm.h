#ifndef COMM_H
#define COMM_H

/* comm
 *
 * $Revision$
 *
 * Author:     R. D. Henderson
 *
 * Description:
 *
 * This is a very basic interface for message passing.  It's meant to be 
 * implemented on top of something else like NX or MPI.  The only purpose
 * for the comm() routines are to simplify the more complicated interface of
 * systems like MPI, and to add a degree of portability.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>

/* General */
int  comm_init  (int *argc, char **argv[]);
int  comm_exit  ();
int  comm_size  ();
int  comm_rank  ();
int  comm_sync  ();

/* Timing routines */
double comm_time ();
double comm_tick ();

/* Testing routines */
int comm_ringtest ();

/* Message Passing (synchronous) */
int  comm_send  (int type, const void *buf, size_t size, int node);
int  comm_recv  (int type, void *buf, size_t size);
int  comm_xchg  (int type, void *sbuf, void *rbuf, size_t size, int node);
int  comm_bcast (void *buf, size_t size, int root);

/* Message Passing (asynchronous) */
int  comm_sendx (int type, const void *buf, size_t size, int node);
int  comm_recvx (int type,       void *buf, size_t size);
int  comm_wait  (int tag);

/* Global reduction */
int  comm_dsum  (int n, const double *u, double *v);
int  comm_dmax  (int n, const double *u, double *v);
int  comm_isum  (int n, const int *u, int *v);
int  comm_imax  (int n, const int *u, int *v);

#endif

