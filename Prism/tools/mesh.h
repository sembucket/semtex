#ifndef MESH_H
#define MESH_H

/* ------------------------------------------------------------------------- 
 * Mesh
 *
 * A "mesh" is the data structure that contains the global quadrature points 
 * for a spectral element grid.  Geometry data is stored in three flat 
 * arrays that have logical dimensions of
 *
 *                 x[K][N][N], y[K][N][N], z[M]
 *
 * where K = number of elements, N = [nr|ns] is the number of polynomial 
 * coefficients per direction in each element, and M = number of physical
 * points in the periodic direction for a 3D computational domain.
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Revision$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "speclib/speclib.h"

typedef struct {
  int    nr;           /* size parameters */
  int    ns;
  int    nz;
  int    nel;
  double *x, *y, *z;   /* storage for the geometry data */
} Mesh;

/* Allocate and load a Mesh */

Mesh* Mesh_alloc (int nr, int ns, int nz, int nel);
void  Mesh_free  (Mesh *m);

Mesh* Mesh_read   (FILE *fp);
int   Mesh_write  (const Mesh *m, FILE *fp);
int   Mesh_interp (Mesh *m, int nr1, int ns1, int nz1, int uniform);
int   Mesh_define (Mesh *m, const Field *u);

#endif
