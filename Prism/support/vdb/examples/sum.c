/*
 * sum test
 *
 * $Id$
 * 
 * Author: R. D. Henderson
 *
 * Read in a spectral element mesh, construct a VDB of mesh points, and
 * then compute the average of function values across element boundaries.
 * In this test an array is created and a value of 1.0 is added for each
 * node in the mesh.  The array is then summed across all shared points and
 * the value at each entry is divided by its multiplicity.  The final 
 * result should be a value of 1.0 stored at each array index.
 *
 * This tests the members() and combine() operations.
 *
 * Copyright (c) 1998-99 R. D. Henderson and Caltech.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "veclib/veclib.h"
#include "comm/comm.h"
#include "../vdb.h"

void average (const vdb_t *vdb, int npts, double *u)
{
  int n;

  /* compute the multiplicity */

  int *mult = (int*) malloc(npts*sizeof(int));
  vdb_members (vdb, mult);

  /* combine data across element boundaries */

  vdb_combine (vdb, u, sizeof(double), 1, comb_dsum);

  /* now divide by the multiplcity */

  for (n = 0; n < npts; n++)
    u[n] *= (1.0/mult[n]);

  free (mult);
}

/* ------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  int nr, ns, nz, nel, iel;
  int npts;
  int nerr;
  int k, p;

  int proc, nprocs;

  double *u;
  double **xyz;
  int    **nodes;

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
      fprintf (stderr, "usage: sum meshfile\n");
    comm_exit();
    exit(-1);
  }

  fp = fopen(argv[1], "r");
  assert(fp);
  fscanf(fp, "%d%d%d%d", &nr, &ns, &nz, &nel);
  fgets (buf, BUFSIZ, fp);

  /* allocate memory for the data, mesh, etc. */

  npts  = nr*ns*nel;
  u     = dvector(0, npts);
  xyz   = dmatrix(0, npts, 0, 2);
  nodes = imatrix(0, nel, 0, nr*ns-1);
  vdb   = vdb_alloc(2, npts, 0, 0);

  /* read the mesh and create the vdb, distributing the elements evenly *
   * over the set of processors.                                        */

  for (k = 0, iel = 0; k < nel; k++) {
    int i;
    if ((k % nprocs) == proc) {
      double coord[2];
      for (i = 0; i < nr*ns; i++) {
	fscanf (fp, "%lf%lf", &coord[0], &coord[1]);
	p = nodes[iel][i] = vdb_djoin(vdb, coord);
	memcpy (xyz[p], coord, sizeof(coord));
      }
      iel++;
    } else
      for (i = 0; i < nr*ns; i++)
	fscanf (fp, "%*lf%*lf");
  }

  npts = vdb_sync(vdb);

  memset (u, 0, npts*sizeof(double));

  for (k = 0; k < iel; k++) {
    for (p = 0; p < nr*ns; p++)
      u[nodes[k][p]] += 1.0;
  }

  average (vdb, npts, u);

  nerr = 0;
  for (p = 0; p < npts; p++)
    if (fabs(u[p]-1.0) > 1.e-6)
      nerr++;

  p = nerr; comm_isum(1, &p, &nerr);

  if (proc == 0) {
    if (nerr==0) 
      fprintf (stderr, "sum_test is OK\n");
    else
      fprintf (stderr, "sum_test failed with %d errors\n", nerr);
  }

  /* Info */

  if (proc == 0) 
    fprintf (stderr, "\npoint distribution:\n");
  for (p = 0; p < nprocs; p++) {
    if (p == proc)
      fprintf (stderr, "\t %2d %4d %6d\n", proc, iel, npts);
    comm_sync();
  }
  if (proc == 0) 
    fprintf (stderr, "\n");

  vdb_stats (vdb, stderr);
  vdb_free  (vdb);

  return comm_exit();
}

  
  
