/*
 * SEM matrix-vector multiply test
 *
 * $Id$
 *
 * Author: R. D. Henderson
 *
 * Read in a spectral element mesh, construct a VDB of mesh points, and
 * then compute the matrix-vector product for the discrete Laplace matrix.
 * This is pretty crude -- in particular, the numerical result is incorrect
 * because none of the geometric scaling is taken care of.  It's just to
 * give some rough timing data.
 *
 * Copyright (c) 1998-99 R. D. Henderson and Caltech.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "cubit/cubit.h"
#include "veclib/veclib.h"
#include "comm/comm.h"
#include "../vdb.h"

static vdb_t *vdb;    /* voxel database */
 
static int nr;        /* mesh resolution (local) */
static int ns;
static int nz;
static int nelmt;     

static int *secr;     /* secretary points */
static int *mult;     /* multiplicity */
static int **key;     /* voxel keys */
  
static double **xyz;  /* mesh */

static int nloop = 100;

static double A_dot (const double *u, const double *v)
{
  const int npts = nr*ns*nelmt;
  double sum = 0.;
  double result;
  int i;

  for (i = 0; i < npts; i++)
    if (secr[i])
      sum += u[i]*v[i];

  comm_dsum (1, &sum, &result);

  return result;
}
  
static void A_mult (const double *u, double *v)
{
  const int nrns = nr*ns;       /* number of points per element */
  const int npts = nr*ns*nelmt; /* number of points on this processor */

  int i, k;
  double **dr, **ds;

  /* Local compute buffers */

  double *uloc = dvector(0, nrns-1);
  double *vloc = dvector(0, nrns-1);
  double *Pr   = dvector(0, nrns-1);
  double *Ps   = dvector(0, nrns-1);

  /* Generate spectral operators */

  getops (nr, 0, 0, &dr, 0);
  getops (ns, 0, 0, &ds, 0);

  /* Clear the result vector */

  dzero (npts, v, 1);

  /* Loop... */

  for (k = 0; k < nelmt; k++) {

    for (i = 0; i < nrns; i++)
      uloc[i] = u[ key[k][i] ];

    dgemm ('T','N', nr, ns, nr, 1., *dr, nr, uloc, nr, 0., Pr, nr); 
    dgemm ('N','N', nr, ns, ns, 1., uloc, nr, *ds, ns, 0., Ps, nr); 

    dgemm ('N', 'N', nr, ns, nr, 1., *dr, nr, Pr, nr, 0., vloc, nr);
    dgemm ('N', 'T', nr, ns, ns, 1., Ps, nr, *ds, ns, 1., vloc, nr);

    for (i = 0; i < nrns; i++)
      v[ key[k][i] ] += vloc[i];
  }

  /* Compute local result vector with remote vectors */

  vdb_combine (vdb, v, sizeof(double), 1, comb_dsum);

  /* Done */

  free_dvector (uloc, 0);
  free_dvector (vloc, 0);
  free_dvector (Pr,   0);
  free_dvector (Ps,   0);

  return;
}

/* ------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  FILE   *fp;
  char   buf[BUFSIZ];
  double *u, *v, result;
  int    npts, i, k, p, proc, nprocs, nglob;

  /* Initialize communications */

  comm_init(&argc, &argv);
  proc   = comm_rank();
  nprocs = comm_size();

  /* Open the data file and read in the points.  Data for each element is *
   * distributed to the processors in order...                            */

  fp = fopen(argv[1], "r");
  assert(fp);
  fscanf(fp, "%d%d%d%d", &nr, &ns, &nz, &nglob);
  fgets (buf, BUFSIZ, fp);

  /* some error checking */

  if (nglob < comm_size()) {
    if (comm_rank()==0) {
      fprintf (stderr, 
	       "mxv: must have at least one element for each processor\n");
    }
    comm_exit();
    exit(0);
  }

  /* allocate memory for the data, mesh, etc. */

  npts = nr*ns*nglob;
  u    = dvector(0, npts-1);
  v    = dvector(0, npts-1);
  xyz  = dmatrix(0, npts-1, 0, 1);

  key  = imatrix(0, nglob, 0, nr*ns-1);
  secr = ivector(0, npts-1);
  mult = ivector(0, npts-1);

  vdb  = vdb_alloc (2, npts, 1 << 10, 1.e-4);

  /* read the mesh and create the vdb, distributing the elements evenly *
   * (but randomly) over the set of processors.                         */

  nelmt = 0;
  for (k = 0; k < nglob; k++) {
    if ((k % nprocs) == proc) {
      double coord[2];
      for (i = 0; i < nr*ns; i++) {
	fscanf (fp, "%lf%lf", &coord[0], &coord[1]);
	key[nelmt][i] = vdb_dadd(vdb, coord);
	memcpy (xyz[key[nelmt][i]], coord, sizeof(coord));
      }
      nelmt++;
    } else
      for (i = 0; i < nr*ns; i++)
	fscanf (fp, "%*lf%*lf");
  }

  npts = vdb_sync(vdb);
  vdb_secr   (vdb,secr);
  vdb_members(vdb,mult);

  /* Generate a solution vector */

  for (i = 0; i < npts; i++)
    u[i] = sin(xyz[i][0]*xyz[i][0] + xyz[i][1]*xyz[i][1]);

  /* CG kernal test */

  A_mult(u,v);
  result = A_dot(u,v);
  if (comm_rank()==0) 
    printf ("result = %g\n", result);

  vdb_stats(vdb,stderr);

  if (comm_rank()==0) 
    fprintf (stderr, "starting timing run:\n");
  
  comm_sync(); {
    double start;
    double stop;
    int n;

    start = dclock();
    for (n = 0; n < nloop; n++) {
      A_mult(u,v);
      result += A_dot(u,v);
    }
    stop = dclock();

    if (comm_rank()==0) 
      fprintf (stderr, "average time = %g sec\n", (stop-start)/nloop);
  }

  return comm_exit();
}

  
  
