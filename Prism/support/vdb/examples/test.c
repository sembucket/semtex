/*
 * vdb test code
 *
 * $Id$
 * 
 * Author: R. D. Henderson
 *
 * This program just reads "commands" from the user and initiates some kind
 * of action in the library.  The commands are:
 *
 *    a <proc> <x> <y> <z>  join the point <x,y,z> to processor <proc>'s VDB
 *    s                     show all vdb's
 *    x                     syncronize all vdb's
 *    q                     quit
 *
 * The file "test.input" provides a sample input file.
 *
 * Copyright (c) 1998-99 R. D. Henderson and Caltech.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "vdb.h"
#include "comm/comm.h"

#define VOXEL_DIM  3
#define VOXEL_MAX  1024

#define DEBUG

typedef enum {
  JOIN,
  DELETE,
  QUERY,
  SHOW,
  SYNC,
  IDLE,
  HELP,
  EXIT
} action_t;

struct {
  char     key;
  action_t action;
} methods[] = {
  { 'j', JOIN },
  { 'd', DELETE },
  { 'q', QUERY },
  { 'h', HELP },
  { 'x', SYNC },
  { 's', SHOW },
  { 'e', EXIT },
};

/* ------------------------------------------------------------------------- */

static void test_master (vdb_t *vdb)
{
  char buf[BUFSIZ];
  int  action = IDLE;
  int  dest;
  double coord[3];

  while (action != EXIT) {
    fputs ("command> ", stdout);
    if (fgets(buf, BUFSIZ, stdin) == NULL)
      sprintf(buf, "q");

    switch (*buf) {
    case 'j':
      action = JOIN;
      sscanf (buf, "%*s%d%lf%lf%lf", &dest, coord, coord+1, coord+2);
      if (dest != 0) {
#ifdef DEBUG
	fprintf (stderr, "sending [%g,%g,%g] -> [%d]\n",
		 coord[0], coord[1], coord[2], dest);
#endif
	comm_send(0, &action, sizeof(action), dest);
	comm_send(1, coord,   sizeof(coord),  dest);
      } else {
	vdb_djoin (vdb, coord);
      }
      break;

    case 'd':
      action = DELETE;
      sscanf (buf, "%*s%d%lf%lf%lf", &dest, coord, coord+1, coord+2);
      if (dest != 0) {
#ifdef DEBUG
	fprintf (stderr, "sending [%g,%g,%g] -> [%d]\n",
		 coord[0], coord[1], coord[2], dest);
#endif
	comm_send(0, &action, sizeof(action), dest);
	comm_send(1, coord,   sizeof(coord),  dest);
      } else {
	vdb_ddelete (vdb, coord);
      }
      break;
      

    case 'q': 
      action = QUERY;
      sscanf (buf, "%*s%d%lf%lf%lf", &dest, coord, coord+1, coord+2);
      if (dest != 0) {
#ifdef DEBUG
	fprintf (stderr, "sending [%g,%g,%g] -> [%d]\n",
		 coord[0], coord[1], coord[2], dest);
#endif
	comm_send(0, &action, sizeof(action), dest);
	comm_send(1, coord,   sizeof(coord),  dest);
      } else {
	printf ("[%d]: [%g,%g,%g] -> %d\n", 
		comm_rank(), coord[0], coord[1], coord[2],
		vdb_dquery(vdb,coord));
      }
      break;
    
    case 's':
      action = SHOW;
      fprintf (stderr, "----- proc 0\n");
      vdb_show (vdb, stderr);
      for (dest = 1; dest < comm_size(); dest++) {
	comm_send (0, &action, sizeof(action), dest);
	comm_recv (0, &action, sizeof(action));
      }
      break;

    case 'x':
      action = SYNC;
      for (dest = 1; dest < comm_size(); dest++)
	comm_send (0, &action, sizeof(action), dest);
      vdb_sync(vdb);
      break;

    case 'e':
      action = EXIT;
      for (dest = 1; dest < comm_size(); dest++) {
	comm_send (0, &action, sizeof(action), dest);
      }
      break;

    default:
      break;
    }
  }

  return;
}
  
/* ------------------------------------------------------------------------- */

static void test_slave (vdb_t *vdb)
{
  int    action = IDLE;
  double coord[3];

  while (action != EXIT) {
    comm_recv (0, &action, sizeof(action));

    switch (action) {
    case JOIN:
      comm_recv (1, coord, sizeof(coord));
      vdb_djoin (vdb, coord);
#ifdef DEBUG
      fprintf (stderr, "[%d]: adding [%g,%g,%g]\n", 
	       comm_rank(), coord[0], coord[1], coord[2]);
#endif
      break;

    case DELETE:
      comm_recv (1, coord, sizeof(coord));
      vdb_ddelete (vdb, coord);
#ifdef DEBUG
      fprintf (stderr, "[%d]: adding [%g,%g,%g]\n", 
	       comm_rank(), coord[0], coord[1], coord[2]);
#endif
      break;

    case QUERY:
      comm_recv (1, coord, sizeof(coord));
      printf ("[%d]: [%g,%g,%g] -> %d\n", 
	      comm_rank(), coord[0], coord[1], coord[2],
	      vdb_dquery(vdb,coord));
      break;
      
    case SHOW:
      fprintf (stderr, "----- proc %d\n", comm_rank());
      vdb_show(vdb, stderr);
      comm_send(0, &action, sizeof(action), 0);
      break;

    case SYNC:
      vdb_sync(vdb);
      break;

    case IDLE:
      break;
    case EXIT:
      break;
    default:
      fprintf (stderr, "[%d]: unknown action\n", comm_rank());
      break;
    }
  }
}
      

/* ------------------------------------------------------------------------- */


main (int argc, char *argv[])
{
  vdb_t *vdb;

  comm_init(&argc, &argv);

  /* allocate the vdb */

  vdb = vdb_alloc(VOXEL_DIM, VOXEL_MAX, 0, FLT_EPSILON);

  if (comm_rank()==0) 
    test_master(vdb);
  else
    test_slave (vdb);

  comm_exit();
}

