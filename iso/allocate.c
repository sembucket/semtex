/*****************************************************************************
 * allocate.c: make storage for user-defined types.
 *
 * : allocate.c,v 2.1 1995/11/08 00:42:04 hmb Exp $
 *****************************************************************************/

#include "iso.h"


int* ivector (int lo,
	      int hi)
/* ------------------------------------------------------------------------- *
 * Allocates an int vector with range [lo..hi].
 * ------------------------------------------------------------------------- */
{
  int* v;
  
  v = (int*) calloc (hi-lo+1, sizeof (int));
  if (!v) message ("ivector", "allocation failure", ERROR);
  return v - lo;
}


complex* cvector (int lo,
		  int hi)
/* ------------------------------------------------------------------------- *
 * Allocates a complex vector with range [lo..hi].
 * ------------------------------------------------------------------------- */
{
  complex* v;
  
  v = (complex*) calloc (hi-lo+1, sizeof (complex));
  if (!v) message ("cvector", "allocation failure", ERROR);
  return v-lo;
}


real* rvector (int lo,
	       int hi)
/* ------------------------------------------------------------------------- *
 * Allocates a real vector with range [lo..hi].
 * ------------------------------------------------------------------------- */
{
  real* v;
  
  v = (real*) calloc (hi-lo+1, sizeof (real));
  if (!v) message ("cvector", "allocation failure", ERROR);
  return v-lo;
}


real* cbox (int ilo,
	    int ihi,
	    int jlo,
	    int jhi,
	    int klo,
	    int khi,
	    CF* c  )
/* ------------------------------------------------------------------------- *
 * Allocates a complex 3-D matrix with ranges [ilo..ihi][jlo..jhi][klo..khi].
 * Returns storage guaranteed to be contiguous.
 * The function returns the start address of the contiguous storage, with
 * the updated argument c a pointer to an array of pointers to an array
 * of pointers to complex.  Returns real*, since we will usually want to
 * use that for fast indexing through every element of storage.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "cbox";
  int       i, j;
  complex  *head;
  real     *fhead;

  /* -- Allocate contiguous storage. */
  head = (complex*)
   calloc ((ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1), sizeof (complex));
  if (!head) message (routine, "contiguous allocation failure", ERROR);
  
  /* -- Allocate pointers to i-direction. */
  *c = (complex ***) calloc (ihi-ilo+1, sizeof (complex **));
  if (!*c) message (routine, "allocation failure 1", ERROR);
  *c -= ilo;

  /* -- Allocate pointers to j-direction. */
  for (i = ilo; i <= ihi; i++) {
    *(*c + i) = (complex**) calloc (jhi-jlo+1, sizeof (complex *));
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


real** cfield (CVF* u)
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

  *u = (CVF) calloc (3, sizeof (CF));
  if (!u) message (routine, "unable to allocate complex vector field", ERROR);
  *u -= 1;   /* indices are [1..3] */

  h = (real**) calloc (3, sizeof (real *));
  if (!h) message (routine, "unable to allocate component handles", ERROR);
  h -= 1;

  h[1] = cbox (0, N-1, 0, N-1, 0, K-1, &((*u)[1]));
  h[2] = cbox (0, N-1, 0, N-1, 0, K-1, &((*u)[2]));  
  h[3] = cbox (0, N-1, 0, N-1, 0, K-1, &((*u)[3]));

  return h;
}


void allocate (CVF*      V    ,
	       CVF*      G    ,
	       const int order,
	       CVF*      WK   ,
	       CF*       F    ,
	       CF*       F_   ,
	       complex** Wtab ,
	       complex** Stab )
/* ------------------------------------------------------------------------- *
 * Use the above routines to get main storage for N-S simulation.
 * ------------------------------------------------------------------------- */
{
  int i;

  cfield (&(*V));
  for (i = 0; i < order; i++) cfield (&G[i]);
  cfield (&(*WK));

  cbox (0, N-1, 0, N-1, 0, K-1, F );
  cbox (0, N-1, 0, N-1, 0, K-1, F_);
  *Wtab = cvector (0, K-1);
  *Stab = cvector (-(N-1), N-1);
}
