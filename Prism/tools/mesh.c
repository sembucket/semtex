/* 
 * Reading / Writing / Processing mesh files
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "speclib/speclib.h"
#include "veclib/veclib.h"
#include "mesh.h"

/* Private Functions */

static void error_msg (char *msg) { 
  fprintf (stderr, "mesh: %s\n", msg); exit(-1); 
}

/* ---------------------------------------------------------------------- */

Mesh* Mesh_alloc (int nr, int ns, int nz, int nel)
{
  Mesh *m = (Mesh*) calloc (1, sizeof(Mesh));
  assert(m);

  m->x = (double*) calloc (nr*ns*nel, sizeof(double));
  m->y = (double*) calloc (nr*ns*nel, sizeof(double));
  m->z = (double*) calloc (nz+1,      sizeof(double));
  assert (m->x && m->y && m->z);

  m->nr  = nr;
  m->ns  = ns;
  m->nz  = nz;
  m->nel = nel;

  return m;
}

void Mesh_free (Mesh *m)
{
  if (m->x) free(m->x);
  if (m->y) free(m->y);
  if (m->z) free(m->z);

  free(m);
}

/* ---------------------------------------------------------------------- *
 * Mesh_read() -- Read a mesh file                                        *
 *                                                                        *
 * This function allocates a mesh structure and loads the mesh from a     *
 * given file.  It returns the new mesh structure.                        *
 * ---------------------------------------------------------------------- */

Mesh *Mesh_read (FILE *fp)
{
  int  nr, ns, nz, nel;
  char buf[BUFSIZ];
  int  i;

  Mesh *m;

  fgets (buf, BUFSIZ, fp);
  if (sscanf (buf, "%d%d%d%d", &nr, &ns, &nz, &nel) != 4)
    error_msg ("unable to read the mesh header");
  
  /* allocate memory */

  m = Mesh_alloc (nr, ns, nz, nel);

  /* read in the mesh (x,y)-first, then (z) */

  for (i = 0; i < nr*ns*nel; i++)
    if (fscanf (fp, "%lf%lf", m->x + i, m->y + i) != 2)
      error_msg ("an error occured while reading the x-y mesh in");
  
  if (nz > 2) {
    for (i = 0; i < nz; i++)
      if (fscanf (fp, "%lf", m->z + i) != 1)
	error_msg ("an error occured while reading the z-mesh in");
  }
  
  return m;
}

/* ------------------------------------------------------------------------- *
 * Mesh_write()                                                              *
 *                                                                           *
 * This function writes a Mesh structure to a mesh file.                     *
 * ------------------------------------------------------------------------- */

int Mesh_write (const Mesh *m, FILE *fp)
{
  const int n  = m->nr * m->ns * m->nel;
  const int nz = m->nz;

  int i;

  fprintf (fp, "%d %d %d %d NR NS NZ NEL\n", 
	   m->nr, m->ns, m->nz, m->nel);

  for (i = 0; i < n; i++)
    if (fprintf (fp, "%#16.10g %#16.10g\n", m->x[i], m->y[i]) < 0)
      return 1;

  if (nz > 1)
    for (i = 0; i <= nz; i++)
      if (fprintf (fp, "%#16.10g\n", m->z[i]) < 0)
	return 1;

  return 0;
}

/* ------------------------------------------------------------------------- *
 * Mesh_define() -- Defines the mesh geometry by copying it from a Field     *
 * ------------------------------------------------------------------------- */

int Mesh_define (Mesh *m, const Field *u)
{
  const int nrns = m->nr * m->ns;
  double    *z   = zmesh(m->nz);

  int k;
  Element *elmt;

  if (m->nr != u->nr || m->ns != u->ns || m->nel != Field_count(u) || 
      m->nz != u->nz)
    error_msg("mesh and field have different sizes");

  for (k = 0, elmt = Field_head(u); elmt != 0; k++, elmt = elmt->next) {
    memcpy (m->x + k*nrns, *elmt->xmesh, nrns*sizeof(double));
    memcpy (m->y + k*nrns, *elmt->ymesh, nrns*sizeof(double));
  }

  memcpy (m->z, z, (m->nz+1)*sizeof(double));
  free   (z);

  return 0;
}
 
/* ---------------------------------------------------------------------- *
 * Mesh_interp() -- Interpolate a mesh to a different NORDER              *
 *                                                                        *
 * This function interpolates a mesh to a finer or coarser resolution.    *
 * If the argument uniform=1 it interpolates to an evenly-spaced set of   *
 * points.                                                                *
 *                                                                        *
 * If any arguments are passed as 0 or the flag UNSET the current values  *
 * are retained.                                                          *
 *                                                                        *
 * NOTE: This function INTERPOLATES the current mesh to a new one, it     *
 *       does not GENERATE a mesh at a different N-order.  This makes a   *
 *       big difference for curved sides.                                 *
 * ---------------------------------------------------------------------- */

static double *do_interp 
  (double **imr, double **itmr, double **ims, double **itms, 
   double *data, int nrp, int nsp, int nr, int ns, int nel)
{
  int     nrns = nr  * ns;
  int     npnp = nrp * nsp;
  int     nxy  = nr  * ns * nel;
  double *new  = dvector (0, npnp * nel - 1), *p = new;
  double *tmp  = dvector (0, nsp*nr);

  int  k;
  for (k = 0; k < nel; k++, data += nrns, p += npnp) {
    mxm (*ims, nsp,  data, ns, tmp, nr );
    mxm ( tmp, nsp, *itmr, nr,   p, nrp);
  }

  free (tmp);
  return new;
}

int Mesh_interp (Mesh *m, int nrp, int nsp, int nzp, int uniform)
{
  int      nr, ns, nel;
  double   **imr, **itmr, **ims, **itms, *zr, *zs, *zrm, *zsm;
  double   *x, *y, *p;

  /* First do the z-mesh */
  
  if (nzp != 0 && (nzp != UNSET)) {
    iparam_set ("NZ", nzp);
    if (m->z) free (m->z);
    m->z = zmesh(m->nz = nzp);
  }

  if (nrp == UNSET && nsp == UNSET) return 0;

  nr     = m->nr;     nrp = nrp ? nrp : nr;
  ns     = m->ns;     nsp = nsp ? nsp : ns;
  nel    = m->nel;
  x      = m->x;
  y      = m->y;

  imr    = dmatrix(0, nrp-1, 0,  nr-1);
  itmr   = dmatrix(0,  nr-1, 0, nrp-1);
  ims    = dmatrix(0, nsp-1, 0,  ns-1);
  itms   = dmatrix(0,  ns-1, 0, nsp-1);

  /* compute the GLL-mesh and the new mesh */

  coef(nr);  getops(nr, &zr, 0, 0, 0);
  coef(ns);  getops(ns, &zs, 0, 0, 0);

  if (option("uniform")) {             /* Interpolate to an even mesh */
    zrm    = dvector(0, nrp-1);
    zrm[0] = -1.; 
    zrm[1] =  2./(nrp-1);
    dramp (nrp, zrm, zrm+1, zrm, 1);     

    zsm    = dvector(0, nsp-1);
    zsm[0] = -1.; 
    zsm[1] =  2./(nsp-1);
    dramp (nsp, zsm, zsm+1, zsm, 1);     
  } else {                             /* Interpolate to another GLL mesh */
    coef (nrp); getops(nrp, &zrm, 0, 0, 0);
    coef (nsp); getops(nsp, &zsm, 0, 0, 0);
  }
  
  /* compute the interpolation matrices */

  igllm (imr, itmr, zr, zrm, nr, nrp);
  igllm (ims, itms, zs, zsm, ns, nsp);

  /* interpolate to the new mesh */

  p    = do_interp (imr, itmr, ims, itms, x, nrp, nsp, nr, ns, nel); free (x); 
  m->x = p;
  p    = do_interp (imr, itmr, ims, itms, y, nrp, nsp, nr, ns, nel); free (y);
  m->y = p;

  m->nr = nrp;
  m->ns = nsp;

  if (option("uniform")) { free (zrm); free (zsm); }

  free_dmatrix (imr , 0, 0);
  free_dmatrix (itmr, 0, 0);
  free_dmatrix (ims , 0, 0);
  free_dmatrix (itms, 0, 0);

  return 0;
}

