/* 
 * hello, world.
 *
 * $Id$
 *
 * Author: R. D. Henderson
 *
 * A test program to make sure that libcomm.a is working
 *
 * Copyright (c) 1998-99 R. D. Henderson and Caltech.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "comm/comm.h"

main (int argc, char *argv[])
{
   int npes, rank, p;

   comm_init(&argc,&argv);

   npes = comm_size();
   rank = comm_rank();

   printf ("Hello from processor %d/%d\n", rank, npes);

   comm_exit();
}

