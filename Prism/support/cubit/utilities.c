/*
 * Cubit Utilities
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * -------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "cubit.h"

#include "veclib/veclib.h"

/* Private functions */

static void 
       edge_init (Edge   *edge, Edge   *link, int id, int nr, int ns),
       vert_init (Vertex *vert, Vertex *link, int id); 

static int     *build_emap    (int nr, int ns);
static double **frame_alloc   (int nr, int ns, int nz, int nel);
static void     frame_project (int nr, int ns, int nz, int nel, double *data,
			       int nx, int ny, double *proj);

/* -------------------------------------------------------------------- * 
 * cubit_err() -- A simpler error message handler                       *
 *                                                                      *
 * This is the error handler used by most of the Cubit routines.  If    *
 * you call it with msg != NULL, it prints your message and dies.  You  *
 * can also write a message into cubit_err_msg[] and call cubit_err     *
 * with a msg = NULL.                                                   *
 * -------------------------------------------------------------------- */

static                                 /* ...... Message Buffers ...... */
char _message_buf[128];                /* error buffer                  */
char *cubit_err_msg = _message_buf;    /* error message                 */
int   cubit_err_num;                   /* error number  (?)             */

void cubit_err (char *msg)
{
  if (msg)
    fprintf (stderr, "cubit: %s\n", msg);
  else
    fprintf (stderr, "cubit: %s\n", cubit_err_msg);

  exit (cubit_err_num);
}

/* Allocate memory for a multi-frame field */

static double **frame_alloc (int nr, int ns, int nz, int nel)
{
  double *dblock, **pblock;
  int i;

  int ntot = nr * ns * nz * nel;  /* Total storage for all frames */

  /* Parallel Transpose:  Make sure the total storage area can be //
  // be partitioned (N N K) (M/P) or (N N K / P) (M).             */

  if (nz > 1) {
    const int nprocs = option("nprocs");
    ntot  = ((ntot + nprocs - 1) / nprocs) * nprocs;
  } 

  pblock = (double**) malloc (ns * nz * nel * sizeof(double*));
  dblock = (double* ) calloc (ntot, sizeof(double));
  assert (pblock != NULL && dblock != NULL);
  
  for (i = 0; i < ns * nel * nz; i++)
    pblock[i] = dblock + i * nr;

  return pblock;
}
  

/* ------------------------------------------------------------------------ *
 * build_emap() - Compute the Elemental Boundary/Interior Scatter Vector    *
 *                                                                          *
 * This function computes the index map which takes a sequentially ordered  *
 * matrix into a boundary/interior ordered matrix.  The map is recomputed   *
 * whenever this function is called with a new value of (nr,ns); otherwise, *
 * the same map computed for the previous (nr,ns) is returned.              *
 * ------------------------------------------------------------------------ */

static int *build_emap (int nr, int ns)
{
  int     pos  = 0;
  static int *emap = NULL, _nr = 0, _ns = 0;
  int i, j, k;

  if (emap == NULL || nr != _nr || ns != _ns) {
    _nr  = nr;
    _ns  = ns;
    emap = (int*) calloc (nr * ns, sizeof(int));
    assert (emap != NULL);

    /* Boundary scatter vector */

    for(j = 0, i = 0; j < nr; j++)
      emap[pos++] = i * nr + j;
    for(j = nr-1, i = 1; i < ns; i++)
      emap[pos++] = i * nr + j;
    for(j = nr-2, i = ns-1; j >= 0; j--)
      emap[pos++] = i * nr + j;
    for(j = 0, i = ns-2; i > 0;i--)
      emap[pos++] = i * nr + j;
    
    /* Interior scatter vector */
    
    for (i = 1;i < ns - 1; i++)
      for (j = 1;j < nr - 1; j++)
	emap[pos++] = i * nr + j;
  }
    
  return emap;
}


/* Project frame data from DATA(nr,ns,nz,K) to PROJ(nx,ny,nz,K) */

static void frame_project
  (int nr, int ns, int nz, int nel, double *data,
   int nx, int ny,                  double *proj)
{
  const int   nrns = nr * ns,
          nxny = nx * ny;
  double *zr, *zs, *zx, *zy;
  double tmp [_MAX_NORDER * _MAX_NORDER];
  int k, m;

  /* Allocate space for the interpolation matrices */

  double **imr   = dmatrix (0, nx-1, 0, nr-1);
  double **itmr  = dmatrix (0, nr-1, 0, nx-1);
  double **ims   = dmatrix (0, ny-1, 0, ns-1);
  double **itms  = dmatrix (0, ns-1, 0, ny-1);

  /* Compute the GLL points */

  coef (nr); getops (nr, &zr, 0, 0, 0);
  coef (ns); getops (ns, &zs, 0, 0, 0);
  coef (nx); getops (nx, &zx, 0, 0, 0);
  coef (ny); getops (ny, &zy, 0, 0, 0);

  /* Compute the interpolation matrices */

  igllm (imr, itmr, zr, zx, nr, nx);
  igllm (ims, itms, zs, zy, ns, ny);

  /* ----- Project ----- */

  for (m = 0; m < nz; m++) {
    for (k = 0; k < nel; k++, data += nrns, proj += nxny) {
      dgemm ('N', 'N', nr, ny, ns, 1., data, nr, *ims, ns, 0., tmp,  nr);
      dgemm ('T', 'N', nx, ny, nr, 1., *imr, nr,  tmp, nr, 0., proj, nx);
    }
  }

  free_dmatrix (imr,  0, 0);
  free_dmatrix (itmr, 0, 0);
  free_dmatrix (ims,  0, 0);
  free_dmatrix (itms, 0, 0);
  return;
}

/* ------------------------------------------------------------------------ *
 * ecopy: special version of "dcopy" with negative skips                    *
 *                                                                          *
 * This does not follow the BLAS conventions on the skip.  A negative skip  *
 * goes backward in memory from the input location ["x" or "y"], not from   *
 * the end of the array ["x" - (n-1)*xskip, "y" - (n-1)*yskip].             *
 *                                                                          *
 * The name "ecopy" implies that this is an edge-of-the-element copy.       *
 * ------------------------------------------------------------------------ */

int ecopy (int n, double *x, int incx, double *y, int incy)
{
  if (incx == 1 && incy == 1)
    memcpy (y, x, n*sizeof(double));
  else 
    while (n--) {
      *y = *x;
      x += incx;
      y += incy;
    }

  return 0;
}

#ifdef DEBUG

/* ------------------------------------------------------------------------ *
 * DEBUGGING ROUTINES                                                       *
 *                                                                          *
 * The following functions display various data structures for debugging    *
 * purposes:                                                                *
 *                                                                          *
 * show_field   - Print an element field                                    *
 * show_matrix  - Print a matrix created with dmatrix                       *
 * show_vector  - Print a vector created with dvector                       *
 * ------------------------------------------------------------------------ */

void show_field (Field *root)
{
  int     nr  = root->nr,
          ns  = root->ns,
          nel = Field_count(root);
  double  **gridx,
          **gridy,
          **field;
  FILE    *fp;
  int i, j, k;

  fp = (option("nprocs") > 1 ? stderr : stdout);

  fprintf(fp, "Showing the field -- %c\n", root->type);
  
  for(k = 0; k < nel; k++) {
    fprintf(fp, "Element = %d\n", root[k].id + 1);

    for(i = ns-1; i >= 0; i--) {
      for(j = 0;j < nr;j++)
	fprintf(fp, "%#7.4lf ", root[k].field[i][j]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  return;
}

void show_matrix(double **matrix,int nrl,int nrh,int ncl,int nch)
{
  int int i, j;
  unsigned int ic = 0,
               nc = 5;

  printf("Showing matrix:\n");
  for(i = nrl; i <= nrh;i++) {
    for(j = ncl, ic = 0; j <= nch; j++, ic++)
      printf("%#10.4lf ", matrix[i][j]);
    putchar('\n');
  }
  return;
}

#define DO_ROW \
   for (i = ncl-1; i <= nch+1; i++)\
     putchar ('-');\
   putchar ('\n');

void show_mbands(double **a, int nrl, int nrh, int ncl, int nch, double eps)
{
  int int i, j;

  DO_ROW;

  for (i = nrl; i <= nrh; i++) {
    printf("|");
    for (j = ncl; j <= nch; j++) {
      if (fabs(a[i][j]) > eps)
	putchar ('x');
      else
	putchar (' ');
    }
    puts ("|");
  }

  DO_ROW;
  return;
}

#undef DO_ROW

void show_vector(double *a, int imin, int imax)
{
  int i;

  printf("Showing vector [%d - %d]:\n", imin, imax);
  for(i = imin; i <= imax; i++)
    printf("%3d: %#14.6g\n", i, a[i]);
  return;
}

void show_ivector(int *a, int imin, int imax)
{
  int i;

  printf("Showing ivector [%d - %d]:\n", imin, imax);
  for(i = imin; i <= imax; i++)
    printf("%3d: %d\n", i, a[i]);
  return;
}

void show_bcs (Bedge Xbc[])
{
  int    i, k, np, skip, nz;
  double *val, *fld;
  Bedge  *bedg;

  printf("Checking boundary system\n");
  
  for(bedg = Xbc; bedg; bedg = bedg->next) {

    printf("Boundary element ID = %d\n",bedg->elmt->id + 1);
    printf("            edge ID = %d\n",bedg->edge->id  + 1);
    printf("Boundary type       = %c\n",bedg->type);
    printf("Boundary info       = ");
    
    if (islower(bedg->type))
      printf("%s\n", bedg->bc.function);
    else
      printf("value\n");

    np    =  bedg->edge->np;
    skip  =  bedg->edge->skip;
    nz    =  bedg->elmt->nz;
    val   =  bedg->bc.value;
    fld   = *bedg->elmt->field + bedg->edge->start;
	
    printf("Boundary values\n");
    for(i = 0;i < np;i++) {
      printf("%2d: ", i);
      for (k = 0; k < nz; k++)
	printf("%#7.4lf ", val[i + k*np]);
      printf("\n");
    }
  }

  return;
}

#endif
