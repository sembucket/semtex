/*****************************************************************************
 * allocate.c: make storage for user-defined types.
 *
 * $Id$
 *****************************************************************************/

#include "globals.h"


int*  ivector (int lo, int hi)
/* ------------------------------------------------------------------------- *
 * Allocates an int vector with range [lo..hi].
 * ------------------------------------------------------------------------- */
{
  int* v;
  
  v = (int*) malloc ((unsigned)(hi-lo+1) * sizeof (int));
  if (!v) message ("ivector", "allocation failure", ERROR);
  return v - lo;
}


complex*  cvector (int lo, int hi)
/* ------------------------------------------------------------------------- *
 * Allocates a complex vector with range [lo..hi].
 * ------------------------------------------------------------------------- */
{
  complex* v;
  
  v = (complex*) malloc((unsigned)(hi-lo+1) * sizeof (complex));
  if (!v) message ("cvector", "allocation failure", ERROR);
  return v-lo;
}


real*  cbox (int ilo,int ihi, int jlo,int jhi, int klo,int khi, CF *c)
/* ------------------------------------------------------------------------- *
 * Allocates a complex 3-D matrix with ranges [ilo..ihi][jlo..jhi][klo..khi].
 * Returns storage guaranteed to be contiguous.
 * The function returns the start address of the contiguous storage, with
 * the updated argument c a pointer to an array of pointers to an array
 * of pointers to complex.  Returns real*, since we will usually want to
 * use that for fast indexing through every element of storage.
 * ------------------------------------------------------------------------- */
{
  char      routine[] = "cbox";
  int       i, j;
  complex  *head;
  real     *fhead;

  /* -- Allocate contiguous storage. */
  head = (complex*)
   malloc((unsigned)(ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1)*sizeof(complex));
  if (!head) message (routine, "contiguous allocation failure", ERROR);
  
  /* -- Allocate pointers to i-direction. */
  *c = (complex ***) malloc((unsigned)(ihi-ilo+1)*sizeof(complex **));
  if (!*c) message (routine, "allocation failure 1", ERROR);
  *c -= ilo;

  /* -- Allocate pointers to j-direction. */
  for (i = ilo; i <= ihi; i++) {
    *(*c + i) = (complex**) malloc ((unsigned)(jhi-jlo+1)*sizeof(complex *));
    if (!(*(*c + i)))  message (routine, "allocation failure 2", ERROR);
    *(*c + i) -= jlo;
  }
  
  /* -- Set pointers to k-direction. */
  for (i = ilo; i <= ihi; i++)
    for (j = jlo; j <= jhi; j++)
      *(*(*c + i) + j) =
	&(head[(i-ilo)*(jhi-jlo+1)*(khi-klo+1) + (j-jlo)*(khi-klo+1)]) - klo;
  
  /* -- Return pointer to start of contiguous storage. */
  fhead = &head[0].Re;
  return fhead;
}


real**  cfield (int* Dim, CVF* u)
/* ------------------------------------------------------------------------- *
 * Allocate a vector[1..3] of cboxes, the indices of which are sized accord-
 * ing to Dim, starting at 0 on each component, i.e., the indices are wave-
 * numbers in the three directions, and the vector indices 1, 2 & 3 are
 * component numbers.  Return a vector of real*, which point in turn at
 * the first storage element of the three cboxes, to be indexed starting at
 * zero.
 * ------------------------------------------------------------------------- */
{
  char    routine[] = "cfield";
  real**  h;

  *u = (CVF) malloc (3 * sizeof (CF));
  if (!u) message (routine, "unable to allocate complex vector field", ERROR);
  *u -= 1;   /* indices are [1..3] */

  h = (real**) malloc (3 * sizeof (real *));
  if (!h) message (routine, "unable to allocate component handles", ERROR);
  h -= 1;

  h[1] = cbox (0, Dim[1]-1, 0, Dim[2]-1, 0, Dim[3]-1, &((*u)[1]));
  h[2] = cbox (0, Dim[1]-1, 0, Dim[2]-1, 0, Dim[3]-1, &((*u)[2]));  
  h[3] = cbox (0, Dim[1]-1, 0, Dim[2]-1, 0, Dim[3]-1, &((*u)[3]));

  return h;
}


void  zero (CVF Z, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Zero all information in Z.
 * ------------------------------------------------------------------------- */
{
  const    int   Npts = 2 * Dim[1] * Dim[2] * Dim[3];
  register int   i;
  register real *h1, *h2, *h3;

  h1 = &Z[1][0][0][0].Re;
  h2 = &Z[2][0][0][0].Re;
  h3 = &Z[3][0][0][0].Re;

  for (i = 0; i < Npts; i++)
    *(h1+i) = (*(h2+i) = (*(h3+i) = 0.0));
}


void allocate_storage (CVF*          V    ,
		       CVF*          G    ,
		       CVF*          G_old,
		       CVF*          WK   ,
		       CF*           F    ,
		       complex**     Wtab ,
		       complex**     Stab ,
		       int*          Dim  )
/* ------------------------------------------------------------------------- *
 * Use the above routines to get main storage for problem.
 * ------------------------------------------------------------------------- */
{
  cfield (Dim, &(*V)    );
  cfield (Dim, &(*G)    );
  cfield (Dim, &(*G_old));
  cfield (Dim, &(*WK)   );

  cbox (0, Dim[1]-1, 0, Dim[2]-1, 0, Dim[3]-1, &(*F));
  *Wtab = cvector (0, Dim[3]-1);
  *Stab = cvector (-(Dim[1]-1), Dim[1]-1);
  
  zero (*G, Dim);
  zero (*G_old, Dim);
}
