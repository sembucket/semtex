/*****************************************************************************
 * family.c: simple storage scheme for families of (0-offset) vectors.
 * Passed a pointer to external storage, return a pointer to internal
 * storage after checking for redundancies.
 *
 * It is *possible* that using families in this way with storage allocated
 * by C++ new can cause problems since technically the outcome of using free
 * on storage allocated by new is undefined.  You can check if this is indeed
 * the cause of some problem by defining DEBUG during compilation, which
 * disables use of families.
 *****************************************************************************/

static char
RCSid[] = "$Id$";

#include <stdio.h>
#include <malloc.h>
#include <alplib.h>

typedef struct ivect {
  int           size ;
  int*          data ;
  int           nrep ;
  struct ivect* next ;
} iVect;

typedef struct dvect {
  int           size ;
  double*       data ;
  int           nrep ;
  struct dvect* next ;
} dVect;

typedef struct svect {
  int           size ;
  float*        data ;
  int           nrep ;
  struct svect* next ;
} sVect;

static iVect* iHead = 0;
static dVect* dHead = 0;
static sVect* sHead = 0;

static int*    iAdopted (const int, const int*);
static double* dAdopted (const int, const double*);
static float*  sAdopted (const int, const float*);

#ifdef DEBUG
  static const int active = 0;	/* -- Disable families. */
#else
  static const int active = 1;
#endif


void iadopt (const int size,
	     int**     vect)
/* ------------------------------------------------------------------------- *
 * If members of *vect have not already been stored then allocate new family
 * storage, load it from *vect, release vect and return pointer to new data
 * area in *vect.  Otherwise return *vect unaltered.  Have also to check that
 * the input pointer is not already in family storage.
 * ------------------------------------------------------------------------- */
{
  iVect* S = 0;
  int*   member;
  
  if (!vect || !*vect) return;

  if (active && (member = iAdopted (size, *vect)) && member != *vect) {
    free (*vect);
    *vect = member;

  } else {
    S = (iVect*) malloc (sizeof (iVect));
    S -> size = size;
    S -> data = *vect;
    S -> nrep = 1;

    if (iHead) S -> next = iHead;
    iHead = S;
  }
}


static int* iAdopted (const int  size,
		      const int* src )
/* ------------------------------------------------------------------------- *
 * Traverse list and see if the entries of vector src (length size) is in
 * storage.  If so, return address of the storage area.  If this is a new
 * occurrence of the entries of src also update nrep, number of replications.
 * ------------------------------------------------------------------------- */
{
  register iVect* p;
  register int    found = 0;

  for (p = iHead; p; p = p -> next) {
    if (p -> size != size)
      continue;
    if (found = src == p -> data)	/* -- Found an alias. */
      break;
    if (found = lisame (size, src, 1, p -> data, 1)) {
      p -> nrep++;
      break;
    }
  }
  
  return found ? p -> data : 0;
}


void iabandon (int** vect)
/* ------------------------------------------------------------------------- *
 * Family deletion operator.  Traverse list looking for a match, decrement
 * nrep and release last copy of internal storage (and list item) if nrep == 0.
 * ------------------------------------------------------------------------- */
{
  register iVect* p;
  register iVect* o = 0;
  register int    found = 0;

  for (p = iHead; p; o = p, p = p -> next)
    if (found = p -> data == *vect) {
      if (--p -> nrep == 0) {
	free (p -> data);
	if   (p == iHead) iHead     = p -> next;
	else              o -> next = p -> next;
	free (p);	
      }
      *vect = 0;
      return;
    }
}


void dadopt (const int size,
	     double**  vect)
/* ------------------------------------------------------------------------- *
 * See comments for iadopt.
 * ------------------------------------------------------------------------- */
{
  dVect*  S = 0;
  double* member;

  if (!vect || !*vect) return;

  if (active && (member = dAdopted (size, *vect)) && member != *vect) {
    free (*vect);
    *vect = member;

  } else {
    S = (dVect*) malloc (sizeof (dVect));
    S -> size = size;
    S -> data = *vect;
    S -> nrep = 1;

    if (dHead) S -> next = dHead;
    dHead = S;
  }
}


static double* dAdopted (const int     size,
			 const double* src )
/* ------------------------------------------------------------------------- *
 * See comments for iAdopted.
 * ------------------------------------------------------------------------- */
{
  register dVect* p;
  register int    found = 0;

  for (p = dHead; p; p = p -> next) {
    if (p -> size != size)
      continue;
    if (found = src == p -> data)	/* -- Found an alias. */
      break;
    if (found = ldsame (size, src, 1, p -> data, 1)) {
      p -> nrep++;
      break;
    }
  }
  
  return found ? p -> data : 0;
}


void dabandon (double** vect)
/* ------------------------------------------------------------------------- *
 * See comments for iabandon.
 * ------------------------------------------------------------------------- */
{
  register dVect* p;
  register dVect* o = 0;
  register int    found = 0;

  for (p = dHead; p; o = p, p = p -> next)
    if (found = p -> data == *vect) {
      if (--p -> nrep == 0) {
	free (p -> data);
	if   (p == dHead) dHead     = p -> next;
	else              o -> next = p -> next;
	free (p);	
      }
      *vect = 0;
      return;
    }
}


void sadopt (const int size,
	     float**   vect)
/* ------------------------------------------------------------------------- *
 * See comments for iadopt.
 * ------------------------------------------------------------------------- */
{
  sVect* S = 0;
  float* member;

  if (!vect || !*vect) return;

  if (active && (member = sAdopted (size, *vect)) && member != *vect) {
    free (*vect);
    *vect = member;

  } else {
    S = (sVect*) malloc (sizeof (sVect));
    S -> size = size;
    S -> data = *vect;
    S -> nrep = 1;

    if (sHead) S -> next = sHead;
    sHead = S;
  }
}


static float* sAdopted (const int    size,
			const float* src )
/* ------------------------------------------------------------------------- *
 * See comments for iAdopted;
 * ------------------------------------------------------------------------- */
{
  register sVect* p;
  register int    found = 0;

  for (p = sHead; p; p = p -> next) {
    if (p -> size != size)
      continue;
    if (found = src == p -> data)	/* -- Found an alias. */
      break;
    if (found = lssame (size, src, 1, p -> data, 1)) {
      p -> nrep++;
      break;
    }
  }
  
  return found ? p -> data : 0;
}


void sabandon (float** vect)
/* ------------------------------------------------------------------------- *
 * See comments for iabandon.
 * ------------------------------------------------------------------------- */
{
  register sVect* p;
  register sVect* o = 0;
  register int    found = 0;

  for (p = sHead; p; o = p, p = p -> next)
    if (found = p -> data == *vect) {
      if (--p -> nrep == 0) {
	free (p -> data);
	if   (p == sHead) sHead     = p -> next;
	else              o -> next = p -> next;
	free (p);	
      }
      *vect = 0;
      return;
    }
}


int FamilySize (int* nint,
		int* ndp ,
		int* nsp )
/* ------------------------------------------------------------------------- *
 * Return total words of family storage, and individual numbers.
 * ------------------------------------------------------------------------- */
{
  int    ni, nd, ns;

  iVect* ip;
  dVect* dp;
  sVect* sp;

  ni = nd = ns = 0;

  for (ip = iHead; ip; ip = ip -> next) ni += ip -> size;
  for (dp = dHead; dp; dp = dp -> next) nd += dp -> size;
  for (sp = sHead; sp; sp = sp -> next) ns += sp -> size;
  
  if (nint) *nint = ni;
  if (ndp ) *ndp  = nd;
  if (nsp ) *nsp  = ns;

  return ni * sizeof (int)    +
         nd * sizeof (double) +
	 ns * sizeof (float);
}



