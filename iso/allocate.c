/*===========================================================================
 * RCS Information:
 * ----------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 *===========================================================================*/


#include "globals.h"

extern void error(string);



ivector ivect(int lo, int hi)
/*===========================================================================*/
/* Allocates an int vector with range [lo..hi].                              */
/*===========================================================================*/
{
  ivector v;
  
  v = (ivector) malloc((unsigned)(hi-lo+1)*sizeof(int));
  if (!v) error("allocation failure in ivector()");
  return v-lo;
}


cvector cvect(int lo, int hi)
/*===========================================================================*/
/* Allocates a complex vector with range [lo..hi].                           */
/*===========================================================================*/
{
  cvector v;
  
  v = (cvector) malloc((unsigned)(hi-lo+1)*sizeof(complex));
  if (!v) error("allocation failure in cvector()");
  return v-lo;
}


float *cbox(int ilo,int ihi, int jlo,int jhi, int klo,int khi, complex_box *c)
/*===========================================================================*/
/* Allocates a complex 3-D matrix with ranges [ilo..ihi][jlo..jhi][klo..khi].*/
/* Returns storage guaranteed to be contiguous.                              */
/* The function returns the start address of the contiguous storage, with    */
/* the updated argument c a pointer to an array of pointers to an array      */
/* of pointers to complex.  Returns float *, since we will usually want to   */
/* use that for fast indexing through every element of storage.              */
/*===========================================================================*/
{
  int      i, j;
  complex  *head;
  float    *fhead;

  /* Allocate contiguous storage */
  head = (complex *)
   malloc((unsigned)(ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1)*sizeof(complex));
  if (!head) error("contiguous allocation failure in cbox");
  
  /* Allocate pointers to i-direction. */
  *c = (complex ***) malloc((unsigned)(ihi-ilo+1)*sizeof(complex **));
  if (!*c) error("allocation failure 1 in cbox()");
  *c -= ilo;

  /* Allocate pointers to j-direction. */
  for (i=ilo; i<=ihi; i++) {
    *(*c + i) = (complex **)malloc((unsigned)(jhi-jlo+1)*sizeof(complex *));
    if (!(*(*c + i))) error("allocation failure 2 in cbox()");
    *(*c + i) -= jlo;
  }
  
  /* Set pointers to k-direction. */
  for (i=ilo; i<=ihi; i++)
    for (j=jlo; j<=jhi; j++)
      *(*(*c + i) + j) =
	&(head[(i-ilo)*(jhi-jlo+1)*(khi-klo+1) + (j-jlo)*(khi-klo+1)])
	  - klo;
  
  /* Return pointer to start of contiguous storage. */
  fhead = &head[0].Re;
  return fhead;
}


component_handle cfield(ivector Dim, complex_vector_field *u)
/*===========================================================================*/
/* Allocate a vector[1..3] of cboxes, the indices of which are sized accord- */
/* ing to Dim, starting at 0 on each component, i.e., the indices are wave-  */
/* numbers in the three directions, and the vector indices 1, 2 & 3 are      */
/* component numbers.  Return a vector of float *, which point in turn at    */
/* the first storage element of the three cboxes, to be indexed starting at  */
/* zero.                                                                     */
/*===========================================================================*/
{
  component_handle     h;

  *u = (complex_vector_field) malloc((unsigned)(3)*sizeof(complex_box));
  if (!u) error("unable to allocate complex vector field");
  *u -= 1;   /* indices are [1..3] */

  h = (float **) malloc(3*sizeof(float *));
  if (!h) error("unable to allocate component handles");
  h -= 1;

  h[1] = cbox(0, Dim[1]-1, 0, Dim[2]-1, 0, Dim[3]-1, &((*u)[1]));
  h[2] = cbox(0, Dim[1]-1, 0, Dim[2]-1, 0, Dim[3]-1, &((*u)[2]));  
  h[3] = cbox(0, Dim[1]-1, 0, Dim[2]-1, 0, Dim[3]-1, &((*u)[3]));

  return h;
}


void zero(/* update */ complex_vector_field Z,
          /* using  */ ivector              Dim)
/*===========================================================================*/
/* Just to ensure that what we get from malloc is zero (I guess calloc would */
/* have done that).                                                          */
/*===========================================================================*/
{
  int   i;
  float *h1, *h2, *h3;

  h1 = &Z[1][0][0][0].Re;
  h2 = &Z[2][0][0][0].Re;
  h3 = &Z[3][0][0][0].Re;

  for (i=0; i<2*Dim[1]*Dim[2]*Dim[3]; i++)
    *(h1+i) = (*(h2+i) = (*(h3+i) = 0.0));
}


void allocate_storage(complex_vector_field *V,
		      complex_vector_field *G,
		      complex_vector_field *G_old,
		      complex_vector_field *WK,
		      complex_box          *F,
		      cvector              *Wtab,
		      cvector              *Stab,
		      ivector               Dim)
/*===========================================================================*/
/* Use the above routines to get main storage for problem.                   */
/*===========================================================================*/
{
  (void) cfield(Dim, &(*V));
  (void) cfield(Dim, &(*G));
  (void) cfield(Dim, &(*G_old));
  (void) cfield(Dim, &(*WK));
  (void) cbox(0, Dim[1]-1, 0, Dim[2]-1, 0, Dim[3]-1, &(*F));
  *Wtab = cvect(0, Dim[3]-1);
  *Stab = cvect(-(Dim[1]-1), Dim[1]-1);
  
  zero(*G, Dim);
  zero(*G_old, Dim);
}
