/*
 * dot-product test
 *
 * $Id$
 *
 * Author: R. D. Henderson
 *
 * Read in a spectral element mesh, construct a VDB of mesh points, and
 * then compute the dot product of a function defined over the mesh.  This
 * is really just a test of computing the secretary points.
 *
 * Copyright (c) 1998-99 R. D. Henderson and Caltech.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "veclib/veclib.h"
#include "comm/comm.h"
#include "../vdb.h"

double dot (const vdb_t *vdb, int npts, double *u, double *v)
{
  double sum = 0.0;
  double result;
  int n;

  /* compute the list of secretary points */

  int *secr = (int*) malloc(npts*sizeof(int));
  vdb_secr (vdb, secr);

  /* compute the local dot product */

  for (n = 0; n < npts; n++)
    if (secr[n])
      sum += u[n]*v[n];
  free (secr);

  /* compute the global sum */

  comm_dsum (1, &sum, &result);

  return result;
}

/* ------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  int nr, ns, nz, nel;
  int npts;
  int k, p;

  int proc, nprocs;

  double sum;
  double *u;
  double *v;
  double **xyz;

  vdb_t  *vdb;

  FILE *fp;
  char buf[BUFSIZ];

  /* Initialize communications */

  comm_init(&argc, &argv);
  proc   = comm_rank();
  nprocs = comm_size();

  /* Open the data file and read in the points.  Data for each element is *
   * distributed to the processors in order...                            */

  if (argc != 2) {
    if (proc == 0)
      fprintf (stderr, "usage: dot_test meshfile\n");
    comm_exit();
    exit(-1);
  }

  fp = fopen(argv[1], "r");
  assert(fp);
  fscanf(fp, "%d%d%d%d", &nr, &ns, &nz, &nel);
  fgets (buf, BUFSIZ, fp);

  /* allocate memory for the data, mesh, etc. */

  npts = nr*ns*nel;
  u    = dvector(0, npts);
  v    = dvector(0, npts);
  xyz  = dmatrix(0, npts, 0, 2);
  vdb  = vdb_alloc(2, npts, 0, 0);

  /* read the mesh and create the vdb, distributing the elements evenly *
   * over the set of processors.                                        */

  for (k = 0; k < nel; k++) {
    int i;
    if ((k % nprocs) == proc) {
      double coord[2];
      for (i = 0; i < nr*ns; i++) {
	fscanf (fp, "%lf%lf", &coord[0], &coord[1]);
	p = vdb_djoin(vdb, coord);
	memcpy (xyz[p], coord, sizeof(coord));
      }
    } else
      for (i = 0; i < nr*ns; i++)
	fscanf (fp, "%*lf%*lf");
  }

  npts = vdb_sync(vdb);

  for (p = 0; p < npts; p++) {
    u[p] = sin(xyz[p][0]*xyz[p][0] + xyz[p][1]*xyz[p][1]);
    v[p] = cos(xyz[p][0]);
  }

  sum = dot(vdb, npts, u, v);
  if (proc == 0) 
    fprintf (stderr, "global sum = %g [expect 1020.68]\n", sum);

  /* Info */

  if (proc == 0) 
    fprintf (stderr, "\npoint distribution:\n");
  for (p = 0; p < nprocs; p++, comm_sync()) {
    if (p == proc)
      fprintf (stderr, "\t %2d %d points\n", proc, npts);
  }
  if (proc == 0) fprintf (stderr, "\n");

  vdb_stats(vdb, stderr);
  vdb_free (vdb);

  return comm_exit();
}

  
  
