/*****************************************************************************
 * mapping.c: build/return integer mapping vectors used to gather/scatter
 * element storage formats from one to another.
 *****************************************************************************/

static char
RCSid_operators[] = "$Id$";

#include <stdio.h>
#include <malloc.h>
#include <alplib.h>
#include <femdef.h>
#include <femlib.h>

typedef struct mapping {
  int             np   ;
  int*            emap ;
  int*            pmap ;
  struct mapping* next ;
} Mapping;

static Mapping* mHead = 0;


void edgemaps (const int nk ,
	       int**     map,
	       int**     inv)
/* ------------------------------------------------------------------------- *
 * An (e)map is an edge-major list of indices, going CCW around the element
 * edges and then traversing the interior in row-major form. 
 * This allows access element storage by traversing edges, tying in with
 * edge-based global numbering.
 *
 * Pmap (the inversion of emap, or inv above) is built afterwards.
 *
 * To convert from row-major storage to boundary-major, then back:
 *   gathr (ntot, value, emap, tmp  );
 *   gathr (ntot, tmp,   pmap, value)  OR  scatr (ntot, tmp, emap, value);
 *
 * If the required maps are "in stock", they are returned from a list,
 * otherwise they are created and added to the list as a first step.
 * NULL input pointers are not assigned a value.
 * ------------------------------------------------------------------------- */
{
  char              routine[] = "edgemaps";
  const int         nk2 = nk * nk;
  register int      found = 0;
  register Mapping* p;

  if (nk < 2)
    message (routine, "input nk < 2", ERROR);


  for (p = mHead; !found && p; p = p -> next) found = nk == p -> np;

  if (!found) {
    register int i, j, k, n;
    const    int nm = nk - 1;
    register int *em, *pm;

    p = (Mapping *) calloc (1, sizeof (Mapping));
    if (mHead) p -> next = mHead;
    mHead = p;

    p -> np   = nk;
    p -> emap = ivector (0, nk2);
    p -> pmap = ivector (0, nk2);

    em = p -> emap;
    pm = p -> pmap;

    /* -- Traverse exterior CCW. */

    k = 0;
    n = 0;
    em[0] = 0;
    for (i = 1; i < nk; i++) em[++k] = n += 1;
    for (i = 1; i < nk; i++) em[++k] = n += nk;
    for (i = 1; i < nk; i++) em[++k] = n -= 1;
    for (i = 1; i < nk; i++) em[++k] = n -= nk;
  
    /* -- Traverse interior in row-major order. */

    for (i = 1; i < nm; i++)
      for (j = 1; j < nm; j++)
	em[k++] = i * nk + j;
  
    /* -- Build pmap. */
  
    for (i = 0; i < nk2; i++) pm[em[i]] = i;
  }

  /* -- p now points to valid storage: return requested operators. */

  if (map) *map = p -> emap;
  if (inv) *inv = p -> pmap;
}
