/*****************************************************************************
 *                      MEMORY ALLOCATION UTILITIES                          *
 *****************************************************************************/

/*------------------*
 * RCS Information: *
 *------------------*/
static char
  RCS_mem[] = "$Id$";


#include <sys/types.h>
#include <malloc.h>
#include <alplib.h>

#define  FREE_ARG   void*





  
complex *cvector(long nl, long nh)
/* ========================================================================= *
 * Allocates a complex vector with subscript range [nl..nh].                 *
 * ========================================================================= */
{
  complex *v;


  v = (complex*) malloc((size_t) ((nh-nl+1)*sizeof(complex)));
  if (!v) message("cvector()", "allocation failure", ERROR);

  return v-nl;
}





double *dvector(long nl, long nh)
/* ========================================================================= *
 * Allocates a double vector with subscript range [nl..nh].                  *
 * ========================================================================= */
{
  double *v;


  v = (double*) malloc((size_t) ((nh-nl+1)*sizeof(double)));
  if (!v) message("dvector()", "allocation failure", ERROR);

  return v-nl;
}





float *svector(long nl, long nh)
/* ========================================================================= *
 * Allocates a float vector with range [nl..nh].                             *
 * ========================================================================= */
{
  float *v;
  

  v = (float*) malloc((size_t) ((nh-nl+1)*sizeof(float)));
  if (!v) message("svector()", "allocation failure", ERROR);

  return v-nl;
}





int *ivector(long nl, long nh)
/* ========================================================================= *
 * Allocates an int vector with subscript range [nl..nh].                    *
 * ========================================================================= */
{
  int *v;

  
  v = (int*) malloc((size_t) ((nh-nl+1)*sizeof(int)));
  if (!v) message("ivector()", "allocation failure", ERROR);

  return v-nl;
}





zomplex *zvector(long nl, long nh)
/* ========================================================================= *
 * Allocates a zomplex vector with subscript range [nl..nh].                 *
 * ========================================================================= */
{
  zomplex *v;


  v = (zomplex*) malloc((size_t) ((nh-nl+1)*sizeof(zomplex)));
  if (!v) message("zvector()", "allocation failure", ERROR);

  return v-nl;
}





complex **cmatrix(long nrl, long nrh, long ncl, long nch)
/* ========================================================================= *
 * Allocate a complex matrix with subscript ranges [nrl..nrh][ncl..nch].     *
 * ========================================================================= */
{
  long      i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  complex **m;


  m = (complex**) malloc((size_t) (nrow*sizeof(complex*)));
  if (!m) message("cmatrix()", "allocation failure 1", ERROR);
  m -= nrl;

  m[nrl] = (complex*) malloc((size_t) (nrow*ncol*sizeof(complex)));
  if (!m[nrl]) message("cmatrix()", "allocation failure 2", ERROR);
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}





double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* ========================================================================= *
 * Allocate a double  matrix with subscript ranges [nrl..nrh][ncl..nch].     *
 * ========================================================================= */
{
  long     i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  double **m;


  m = (double**) malloc((size_t) (nrow*sizeof(double*)));
  if (!m) message("dmatrix()", "allocation failure 1", ERROR);
  m -= nrl;

  m[nrl] = (double*) malloc((size_t) (nrow*ncol*sizeof(double)));
  if (!m[nrl]) message("dmatrix()", "allocation failure 2", ERROR);
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}





float **smatrix(long nrl, long nrh, long ncl, long nch)
/* ========================================================================= *
 * Allocate a float matrix with subscript ranges [nrl..nrh][ncl..nch].       *
 * ========================================================================= */
{
  long    i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  float **m;


  m = (float**) malloc((size_t) (nrow*sizeof(float*)));
  if (!m) message("smatrix()", "allocation failure 1", ERROR);
  m -= nrl;

  m[nrl] = (float*) malloc((size_t) (nrow*ncol*sizeof(float)));
  if (!m[nrl]) message("smatrix()", "allocation failure 2", ERROR);
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}





int **imatrix(long nrl, long nrh, long ncl, long nch)
/* ========================================================================= *
 * Allocate an int matrix with subscript ranges [nrl..nrh][ncl..nch].        *
 * ========================================================================= */
{
  long    i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  int   **m;


  m = (int**) malloc((size_t) (nrow*sizeof(int*)));
  if (!m) message("imatrix()", "allocation failure 1", ERROR);
  m -= nrl;

  m[nrl] = (int*) malloc((size_t) (nrow*ncol*sizeof(int)));
  if (!m[nrl]) message("imatrix()", "allocation failure 2", ERROR);
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}





zomplex **zmatrix(long nrl, long nrh, long ncl, long nch)
/* ========================================================================= *
 * Allocate a zomplex matrix with subscript ranges [nrl..nrh][ncl..nch].     *
 * ========================================================================= */
{
  long      i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  zomplex **m;


  m = (zomplex**) malloc((size_t) (nrow*sizeof(zomplex*)));
  if (!m) message("cmatrix()", "allocation failure 1", ERROR);
  m -= nrl;

  m[nrl] = (zomplex*) malloc((size_t) (nrow*ncol*sizeof(zomplex)));
  if (!m[nrl]) message("cmatrix()", "allocation failure 2", ERROR);
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}





complex ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* ========================================================================= *
 * Allocate a complex 3-tensor with ranges [nrl..nrh][ncl..nch][ndl..ndh].   *
 * ========================================================================= */
{
  long       i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  complex ***t;


  t = (complex***) malloc((size_t) (nrow*sizeof(complex**)));
  if (!t) message("s3tensor()", "allocation failure 1", ERROR);
  t -= nrl;

  t[nrl] = (complex**) malloc((size_t) (nrow*ncol*sizeof(complex*)));
  if (!t[nrl]) message("s3tensor()", "allocation failure 2", ERROR);
  t[nrl] -= ncl;

  t[nrl][ncl] = (complex*) malloc((size_t) (nrow*ncol*ndep*sizeof(complex)));
  if (!t[nrl][ncl]) message("s3tensor()", "allocation failure 3", ERROR);
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}





double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* ========================================================================= *
 * Allocate a double 3-tensor with ranges [nrl..nrh][ncl..nch][ndl..ndh].    *
 * ========================================================================= */
{
  long      i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  double ***t;


  t = (double***) malloc((size_t) (nrow*sizeof(double**)));
  if (!t) message("s3tensor()", "allocation failure 1", ERROR);
  t -= nrl;

  t[nrl] = (double**) malloc((size_t) (nrow*ncol*sizeof(double*)));
  if (!t[nrl]) message("s3tensor()", "allocation failure 2", ERROR);
  t[nrl] -= ncl;

  t[nrl][ncl] = (double*) malloc((size_t) (nrow*ncol*ndep*sizeof(double)));
  if (!t[nrl][ncl]) message("s3tensor()", "allocation failure 3", ERROR);
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}





float ***s3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* ========================================================================= *
 * Allocate a float 3-tensor with ranges [nrl..nrh][ncl..nch][ndl..ndh].     *
 * ========================================================================= */
{
  long     i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  float ***t;


  t = (float***) malloc((size_t) (nrow*sizeof(float**)));
  if (!t) message("s3tensor()", "allocation failure 1", ERROR);
  t -= nrl;

  t[nrl] = (float**) malloc((size_t) (nrow*ncol*sizeof(float*)));
  if (!t[nrl]) message("s3tensor()", "allocation failure 2", ERROR);
  t[nrl] -= ncl;

  t[nrl][ncl] = (float*) malloc((size_t) (nrow*ncol*ndep*sizeof(float)));
  if (!t[nrl][ncl]) message("s3tensor()", "allocation failure 3", ERROR);
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}





int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* ========================================================================= *
 * Allocate an int 3-tensor with ranges [nrl..nrh][ncl..nch][ndl..ndh].      *
 * ========================================================================= */
{
  int     i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  int  ***t;


  t = (int***) malloc((size_t) (nrow*sizeof(int**)));
  if (!t) message("s3tensor()", "allocation failure 1", ERROR);
  t -= nrl;

  t[nrl] = (int**) malloc((size_t) (nrow*ncol*sizeof(int*)));
  if (!t[nrl]) message("s3tensor()", "allocation failure 2", ERROR);
  t[nrl] -= ncl;

  t[nrl][ncl] = (int*) malloc((size_t) (nrow*ncol*ndep*sizeof(int)));
  if (!t[nrl][ncl]) message("s3tensor()", "allocation failure 3", ERROR);
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}





zomplex ***z3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* ========================================================================= *
 * Allocate a zomplex 3-tensor with ranges [nrl..nrh][ncl..nch][ndl..ndh].   *
 * ========================================================================= */
{
  long       i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  zomplex ***t;


  t = (zomplex***) malloc((size_t) (nrow*sizeof(zomplex**)));
  if (!t) message("s3tensor()", "allocation failure 1", ERROR);
  t -= nrl;

  t[nrl] = (zomplex**) malloc((size_t) (nrow*ncol*sizeof(zomplex*)));
  if (!t[nrl]) message("s3tensor()", "allocation failure 2", ERROR);
  t[nrl] -= ncl;

  t[nrl][ncl] = (zomplex*) malloc((size_t) (nrow*ncol*ndep*sizeof(zomplex)));
  if (!t[nrl][ncl]) message("s3tensor()", "allocation failure 3", ERROR);
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}





void free_cvector(complex *v, long nl)
/* ========================================================================= *
 * Frees a complex vector allocated by cvector().                            *
 * ========================================================================= */
{
  free((FREE_ARG) (v+nl));
}





void free_dvector(double *v, long nl)
/* ========================================================================= *
 * Frees a double vector allocated by dvector().                             *
 * ========================================================================= */
{
  free((FREE_ARG) (v+nl));
}




void free_svector(float *v, long nl)
/* ========================================================================= *
 * Frees a float vector allocated by svector().                              *
 * ========================================================================= */
{
  free((FREE_ARG) (v+nl));
}





void free_ivector(int *v, long nl)
/* ========================================================================= *
 * Frees an int vector allocated by ivector().                               *
 * ========================================================================= */
{
  free((FREE_ARG) (v+nl));
}





void free_zvector(zomplex *v, long nl)
/* ========================================================================= *
 * Frees a complex vector allocated by zvector().                            *
 * ========================================================================= */
{
  free((FREE_ARG) (v+nl));
}





void free_cmatrix(complex **m, long nrl, long ncl)
/* ========================================================================= *
 * Frees a complex matrix allocated with cmatrix().                          *
 * ========================================================================= */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}





void free_dmatrix(double **m, long nrl, long ncl)
/* ========================================================================= *
 * Frees a double matrix allocated with dmatrix().                           *
 * ========================================================================= */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}





void free_smatrix(float **m, long nrl, long ncl)
/* ========================================================================= *
 * Frees a float matrix allocated with smatrix().                            *
 * ========================================================================= */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}





void free_imatrix(int **m, long nrl, long ncl)
/* ========================================================================= *
 * Frees an int matrix allocated with imatrix().                             *
 * ========================================================================= */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}





void free_zmatrix(zomplex **m, long nrl, long ncl)
/* ========================================================================= *
 * Frees a zomplex matrix allocated with zmatrix().                          *
 * ========================================================================= */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}





void free_c3tensor(complex ***t, long nrl, long ncl, long ndl)
/* ========================================================================= *
 * Frees a complex 3-tensor allocated with c3tensor().                       *
 * ========================================================================= */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}





void free_d3tensor(double ***t, long nrl, long ncl, long ndl)
/* ========================================================================= *
 * Frees a double 3-tensor allocated with d3tensor().                        *
 * ========================================================================= */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}





void free_s3tensor(float ***t, long nrl, long ncl, long ndl)
/* ========================================================================= *
 * Frees a float 3-tensor allocated with s3tensor().                         *
 * ========================================================================= */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}





void free_i3tensor(int ***t, long nrl, long ncl, long ndl)
/* ========================================================================= *
 * Frees an int 3-tensor allocated with i3tensor().                          *
 * ========================================================================= */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}





void free_z3tensor(zomplex ***t, long nrl, long ncl, long ndl)
/* ========================================================================= *
 * Frees a zomplex 3-tensor allocated with z3tensor().                       *
 * ========================================================================= */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}
