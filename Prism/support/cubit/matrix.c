/*
 * Spectral Element Stiffness Matrix
 *
 * 
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * -------------------------------------------------------------------- */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "vdb/vdb.h"
#include "veclib/veclib.h"

#include "cubit.h"
#include "vertex.h"
#include "edge.h"
#include "alias.h"
#include "element.h"

/* Private functions */

static void build_CG    (Matrix *A, Mesh *mesh);
static void build_map   (Matrix *A, Mesh *mesh);
static void build_solve (Matrix *A, Mesh *mesh);
static void build_mass  (Matrix *A, Mesh *mesh);
static void build_secr  (Matrix *A, Mesh *mesh);
static void build_diag  (Matrix *A, Mesh *mesh);

/* ------------------------------------------------------------------------- */

Matrix *Matrix_alloc (void)
{
  Matrix *A;

  A = (Matrix*) calloc (1, sizeof(Matrix));
  assert(A);

  return A;
}

void Matrix_free (Matrix *A) 
{
  /* free VDBs */
  /* free matrix storage */
  /* free preconditioner storage */
  
  free(A);
}

/* ------------------------------------------------------------------------- *
 * Matrix_update() -- Update the stiffness matrix                            *
 *                                                                           *
 * After updating the mesh (or creating one) the components of the stiffness *
 * matrix need to be updated.  This function updates the various size param- *
 * eters, computes the global mass matrix, and computes any undefined matrix *
 * elements.                                                                 *
 * ------------------------------------------------------------------------- */

void Matrix_update (Matrix *A, Mesh *mesh)
{
  A->nvert     = vdbnmember (mesh->vertx, NULL);
  A->nedge     = vdbnmember (mesh->edges, NULL);

  A->bpts      = A->nvert + A->nedge * (mesh->nr - 2);
  A->ipts      = mesh->elements * SQR(mesh->nr - 2);

  A->singular  = UNSET;
  A->constant  = dparam("LAMBDA");
  A->condition = 0.;

  build_map   (A, mesh);
  build_solve (A, mesh);
  build_secr  (A, mesh);
  build_mass  (A, mesh);
  build_CG    (A, mesh);
  build_diag  (A, mesh);

  return;
}

/* ------------------------------------------------------------------------- */

static void build_secr (Matrix *A, Mesh *mesh)
{
  const int nvert = A->nvert;
  const int nedge = A->nedge;
  int*  secr;

  if (A->secr) 
    free (A->secr);
  A->secr = secr = (int*) calloc (nvert + nedge, sizeof(int));
  assert (secr);

  vdbsecr (mesh->vertx, secr);
  vdbsecr (mesh->edges, secr + nvert);
  return;
}

/* ------------------------------------------------------------------------- */

static void build_solve (Matrix *A, Mesh *mesh)
{
  const int nvert = A->nvert;
  const int nedge = A->nedge;
  const int bpts  = nvert + nedge;
  int i;

  int     *solve;
  Element *elmt;
  Edge    *edge;
  Vertex  *vert;

  if (A->solve) 
    free (A->solve);
  A->solve = solve = (int*) malloc (bpts*sizeof(int));

  /* Turn solve flag "on" for all points initially */

  for (i = 0; i < bpts; i++)
    solve[i] = 1;
 
  /* Now check edge by edge */

  for (elmt = mesh->head; elmt; elmt = elmt->next) {

    /* Aliased points are NOT solved for */

    for (i = 0, vert = elmt->vert_list; i < 4; i++) 
      if (vert[i].alias)
	solve[vert[i].key] = 0;

    for (i = 0, edge = elmt->edge_list; i < 4; i++) {
      if (edge[i].bc) {
	Vertex *a  = edge[i].a;
	Vertex *b  = edge[i].b;
	BC     *bc = edge[i].bc;
	assert (a->alias == NULL && b->alias == NULL);

	/* Take the logcal AND with return value of (user_bc), which *
	 * must return 0 (ESSENTIAL) or 1 (NATURAL) depending on the *
	 * type of boundary condition.                               */

	if (user_bc != NULL) {
	  solve[a->key] &= (*user_bc)(bc, elmt, a->id, -1, NULL);
	  solve[b->key] &= (*user_bc)(bc, elmt, b->id, -1, NULL);
	  solve[edge[i].key + nvert] &= 
	    (*user_bc)(bc, elmt, -1, edge[i].id, NULL);
	} else {
	  solve[a->key] = 0;
	  solve[b->key] = 0;
	  solve[edge[i].key + nvert] = 0;
	}
      }
    }
  }

  vdbcombine 
    (mesh->vertx, solve,         vdbiand, sizeof(int), 1, sizeof(int));
  vdbcombine 
    (mesh->edges, solve + nvert, vdbiand, sizeof(int), 1, sizeof(int));

  return;
}

/* ------------------------------------------------------------------------- */

static void build_map (Matrix *A, Mesh *mesh)
{
  Element *elmt;
  Edge    *edge;
  Vertex  *vert;

  int np    = mesh->nr;
  int nb    = np * 4 - 4;
  int n     = np - 1;
  int nvert = A->nvert;
  int i;

  if (A->map) 
    free (A->map);
  A->map = (struct map*) calloc (_MAX_NEL, sizeof(struct map));

  for (elmt = mesh->head; elmt; elmt = elmt->next) {

    int*  npts  = A->map[elmt->id].npts  
                = (int*) calloc (nb, sizeof(int));
    int** index = A->map[elmt->id].index 
                = (int**) calloc (nb, sizeof(int*));
    double **Z  = A->map[elmt->id].Z     
                = (double**) calloc (nb, sizeof(double*));

    for (vert = elmt->vert_list; vert; vert = vert->next) {
      const int ip = vert->id * (np - 1);
      
      if (vert->alias == NULL) {
	npts [ip]    = 1;
	index[ip]    = (int*) calloc (1, sizeof(int));
	index[ip][0] = vert->key;
      } else {
	npts [ip]    = np;
	index[ip]    = (int*) calloc (np, sizeof(int));

	index[ip][0] = vert->alias->key.a;
	index[ip][n] = vert->alias->key.b;
	Z    [ip]    = vert->alias->type == 1 ? vert->alias->Z [n] : 
	                                        vert->alias->Z [0] ;
	for (i = 1; i < n; i++)
	  index[ip][i] = vert->alias->key.edge * (np-2) + nvert + i-1;
      }
    }

    for (edge = elmt->edge_list; edge; edge = edge->next) {
      const int dir = edge->dir;
      int   ip      = edge->id  * (np - 1) + 1;
      int   istart  = edge->key * (np - 2) + nvert;

      if (dir < 0) istart += (np-3);

      if (edge->alias == NULL || edge->alias->type == 0) {
	for (i = 0; i < np-2; i++, ip++) {
	  index[ip]    = (int*) calloc (1, sizeof(int));
	  index[ip][0] = istart + i*dir;
	  npts [ip]    = 1;
	}
      } else {
	index[ip]    = (int*) calloc (np, sizeof(int));
	index[ip][0] = edge->alias->key.a;
	index[ip][n] = edge->alias->key.b;
	for (i = 1; i < n; i++)
	  index[ip][i] = edge->alias->key.edge * (np-2) + nvert + i-1;

	for (i = 0; i < np-2; i++) {
	  npts [ip + i] = np;
	  index[ip + i] = index[ip];
	  Z    [ip + i] = edge->alias->Z[n-i-1];
	}
      }
    }
  }

  return;
}

static void build_mass (Matrix *A, Mesh *mesh)
{
  int bpts = A->bpts;

  const int np    = mesh->head->nr-2;
  const int nb    = 4 * np + 4;
  const int nvert = A->nvert;

  double  *mass, *minv, mbuf[_MAX_NB];
  Element *elmt;
  int i, j, k;

  if (A->mass) free (A->mass);

  A->mass    = mass = dvector (0, 2*bpts);
  A->massinv = minv = mass + bpts;

  memset (mass, '\0', bpts*sizeof(double));

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    int    *npts   = A->map[elmt->id].npts;
    int    **index = A->map[elmt->id].index;
    double **Z     = A->map[elmt->id].Z;
    
    for (i = 0; i < nb; i++)
      mbuf[i] = (*elmt->mass)[elmt->emap[i]];
    
    for (i = 0; i < nb; i++) {
      if (npts[i] == 1) {
	assert(index[i][0] < bpts);
	mass[index[i][0]] += mbuf[i];

	/* Compute the product:  M = trans(Z) B Z */

      } else {
	for (j = 0; j < npts[i]; j++) {
	  for (k = 0; k < npts[i]; k++)
	    mass[index[i][j]] += Z[i][j] * mbuf[i] * Z[i][k];
	}
      }
    }
  }
  
  /* Sum across the processors */

  vdbcombine 
    (mesh->vertx, mass, vdbdsum, sizeof(double), 1, sizeof(double));
  vdbcombine 
    (mesh->edges, mass+nvert, vdbdsum, sizeof(double), np, np*sizeof(double));
  
  /* Invert it */

  for (i = 0; i < bpts; i++) {
    if (mass[i] > 0.)    
      minv[i] = 1./mass[i];
  }

  return;
}

/* ------------------------------------------------------------------------- */

static void build_CG (Matrix *A, Mesh *mesh)
{
  const int  np   = mesh->head->nr;
  const int  nrns = np*np;
  Matrix_CG *CG   = A->CG;
  Element   *elmt;
  int   i;
  
  if (!CG) CG = (Matrix_CG*) calloc (_MAX_NEL, sizeof(Matrix_CG));

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    const int k = elmt->id;

    /* Check to see if this matrix set can be skipped */

    if (CG[k].family == elmt->family && 
	CG[k].key    == elmt->key    && CG[k].g1 != NULL) 
      continue;
  
    /* Allocate and compute the element stiffness matrix elements */
    
    CG[k].g1 = dmatrix (0, np-1, 0, np-1);   /* Allocate Memory */
    CG[k].g2 = dmatrix (0, np-1, 0, np-1);
    CG[k].g3 = dmatrix (0, np-1, 0, np-1);
    
    CG[k].family = elmt->family;
    CG[k].key    = elmt->key;
    
    if (elmt->xs)
      for (i = 0; i < nrns; i++) {
	(*CG[k].g1)[i] = SQR((*elmt->xs)[i]) + SQR((*elmt->ys)[i]);
	(*CG[k].g3)[i] = (*elmt->xr)[i] * (*elmt->xs)[i];
      }
    else
      for (i = 0; i < nrns; i++) {
	(*CG[k].g1)[i] = SQR((*elmt->ys)[i]);
	(*CG[k].g3)[i] = 0.;
      }
    
    if (elmt->yr)
      for (i = 0; i < nrns; i++) {
	(*CG[k].g2)[i]  = SQR((*elmt->xr)[i]) + SQR((*elmt->yr)[i]);
	(*CG[k].g3)[i] += (*elmt->yr)[i] * (*elmt->ys)[i];
      }
    else
      for (i = 0; i < nrns; i++)
	(*CG[k].g2)[i]  = SQR((*elmt->xr)[i]);
    
    for (i = 0; i < nrns; i++) {
      const double scale = (*elmt->mass)[i] / SQR((*elmt->jac)[i]);
      
      (*CG[k].g1)[i] *=  scale;
      (*CG[k].g2)[i] *=  scale;
      (*CG[k].g3)[i] *= -scale;
    }
    
    if (!(elmt->xs) && !(elmt->yr)) { 
      free_dmatrix (CG[k].g3, 0, 0); CG[k].g3 = NULL; 
    }
  }
  
  A->CG = CG;
}

static void build_diag (Matrix *A, Mesh *mesh)
{
  const int nr    = mesh->nr;
  const int ns    = mesh->ns;
  const int nrns  = mesh->nr*mesh->ns;
  const int nb    = (nr + ns - 2) << 1;
  const int ni    = nrns - nb;
  const int npts  = A->bpts + A->ipts;
  double lambda   = A->constant;

  double **g1, **g2, **g3, **drt, **dst;
  int    *emap, *bmap;
  int i, j, k, id;

  double *M, *dM, *Mi;
  double Mk[_MAX_NORDER*_MAX_NORDER], tmp[_MAX_NORDER], mbuf[_MAX_NB];

  Element *elmt;

  getops (nr, 0, 0, 0, &drt);
  getops (ns, 0, 0, 0, &dst);

  if (A->diag) free (A->diag);

  M  = A->diag = (double*) malloc(npts*sizeof(double));
  Mi = M + A->bpts;
  memset (M, 0, npts*sizeof(double));

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    int    *npts   = A->map[elmt->id].npts;
    int    **index = A->map[elmt->id].index;
    double **Z     = A->map[elmt->id].Z;

    id = elmt->id;
    g1 = A->CG[id].g1;
    g2 = A->CG[id].g2;
    g3 = A->CG[id].g3;
    memset (dM = Mk, 0, nrns*sizeof(double));

    for (i = 0; i < ns; i++) {
      for (j = 0; j < nr; j++, dM++) {
	dvmul(nr, drt[j], 1, drt[j], 1, tmp, 1);
	*dM  = ddot(nr,  g1 [i], 1, tmp, 1);	  
	dvmul(ns, dst[i], 1, dst[i], 1, tmp, 1);
	*dM += ddot(ns,  *g2+j, nr, tmp, 1);

	if (g3)
	  *dM += 2. * g3[i][j] * drt[j][j] * dst[i][i];
      }
    }

    for (i = 0; i < nb; i++)
      mbuf[i] = Mk[elmt->emap[i]];
    for (i = 0; i < nb; i++) {
      if (npts[i] == 1) {
	assert(index[i][0] < A->bpts);
	M[index[i][0]] += mbuf[i];
      } else {
	for (j = 0; j < npts[i]; j++) {
	  for (k = 0; k < npts[i]; k++)
	    M[index[i][j]] += Z[i][j] * mbuf[i] * Z[i][k];
	}
      }
    }
    

    /* NB: This has to be Mi + ni*elmt->id because the id's are not     *
     * contiguous...it depends to the currently active set of elements. *
     * There will be zero blocks in the preconditioner that never get   *
     * accessed.                                                        */

    dgathr (ni, Mk, elmt->emap + nb, Mi + ni*elmt->id);
  }

  /* Sum across the processors */

#if 0
  vdbcombine
    (mesh->vertx, M, vdbdsum, sizeof(double), 1, sizeof(double));
  vdbcombine
    (mesh->edges, M+A->nvert, vdbdsum, sizeof(double), np, np*sizeof(double));
#endif

  /* Invert and make sure it's postive definite... */

  for (i = 0; i < npts; i++)
    M[i] = fabs(M[i]);
  
  for (i = 0; i < A->bpts; i++)
    if (M[i] != 0.)
      M[i] = 1./M[i];

  for (elmt = mesh->head; elmt; elmt = elmt->next)
    dvrecp (ni, Mi + ni*elmt->id, 1, Mi + ni*elmt->id, 1);

#if 0
  printf ("Preconditioner:\n");
  for (i = 0; i < npts; i++)
    printf ("M[%3d] = %g\n", i, M[i]);
#endif
}
