/* 
 * SPECTRAL ELEMENT FAMILY SYSTEM                       
 *                                                                          
 * The following functions implement a family system for a spectral element 
 * mesh.  Any elements which share the same size and shape also share the   
 * same geometric factors and coefficient matrix.  To put it another way,  
 * any two elements with the same Jacobian for the transformation from 
 * physical to computational space are identical in the new space.  Thus, by 
 * keeping up with the elements being added to the grid, it is possible to 
 * significantly reduce the total amount of storage by simply giving these 
 * elements pointers to a single area of memory where the geometric factors 
 * are stored.  The following functions provide a way to manipulate the 
 * family system without any "outside" knowledge of how it is constructed.
 *
 * Here is a list of the Family_* functions:
 *
 *     get             - Get the family for a particular element
 *     create          - Create a new family
 *     add             - Add an element to the given family         
 *
 *     set             - Mark a family as having been processed
 *     reset           - Reset all families                   
 *     count           - Count the number of element families in a Element       
 *                                                                          
 * The functions set() and reset() allow you to perform a "family loop."  
 * The lists and effects of families are otherwise invisible outside of this 
 * file.
 *                                                                          
 * A call to the function Family_disable() turns off the use of families.
 *                                                                          
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "cubit.h"

#include "veclib/veclib.h"

#define GEOMEMBERS  10    /* # of geometry factors in a family */

/* Externals */

static int     use_families = 1;
static Family *family_list  = (Family*) NULL;

/* Private Functions */

static int  inFamily  (Element *, Family *);

/* ------------------------------------------------------------------------ */

Family *Family_create (Element *U)
{
  Family *fam  = (Family *) calloc (1, sizeof(Family));
  
  /* Set family characteristics */
  
  fam->nr          = U->nr;
  fam->ns          = U->ns;
  fam->members     = 1; 
  fam->set         = 0;

  /* Save the geometry data for this family */

  fam->xr   = U->xr;
  fam->xs   = U->xs;
  fam->yr   = U->yr;
  fam->ys   = U->ys;
  fam->jac  = U->jac;
  fam->rx   = U->rx;
  fam->ry   = U->ry;
  fam->sx   = U->sx;
  fam->sy   = U->sy;
  fam->mass = U->mass;

  /* Add to the family list */

  fam->id          = (family_list) ? family_list->id + 1 : 0;
  fam->parent      = U;
  fam->next        = family_list;  
  family_list      = fam;

  if (!use_families) fam->id = U->id;

  return fam;
}


/* The following function should be called once before any  *
 * pre-processing to disable the use of element families.   */

void Family_disable() { use_families = 0; return; }


/* This just adds the given element to an existing family.  Returns *
 * the new number of members in the family.                         */

int Family_add (Family *fam, Element *U)
{
  U->xr   = fam->xr;
  U->xs   = fam->xs;
  U->yr   = fam->yr;
  U->ys   = fam->ys;
  U->jac  = fam->jac;
  U->rx   = fam->rx;
  U->ry   = fam->ry;
  U->sx   = fam->sx;
  U->sy   = fam->sy;
  U->mass = fam->mass;

  return ++(fam->members);
}

/* Returns a pointer to the family for a given element.  If the family    *
 * system is disabled, this function just matches ID numbers.  Otherwise  *
 * it matches the perimeter of the element and the family.                *
 *                                                                        *
 * Returns NULL of no family is found.                                    */

Family *Family_get (Element *U)
{
  Family *fp;
  for (fp = family_list; fp; fp = fp->next)
    if (inFamily (U, fp)) break;
  return fp;
}

void Family_reset (void)
{
  Family *fam;
  for (fam = family_list; fam; fam = fam->next) fam->set = 0;
  return;
}


/* Mark this family as set...for doing loops over families. *
 * Returns the number of times the family has been set.     */

int Family_set (Family *fam) { return ++(fam->set); }

/* This returns the number of families in a linked element list. *
 * Can be less than the number of members of the family list.    */

int Family_count (Element *U)
{
  int      nel   = 0;
  int      count = 0;
  Element *elmt;
  Family  *fp;
  int k;

  for (elmt = U; elmt; elmt = elmt->next)
    nel++;

  if (!use_families) return nel;

  Family_reset();
  for (k = 0; k < nel; k++) {
    if ((fp = Family_get(& U[k])) && !fp->set) {
      count++;
      Family_set (fp);
    }
  }

  return count;
}

#ifdef DEBUG

int Family_check (Element *U)
{
  int  nel   = Field_count(U);
  int *klist = ivector (0, nel);
  Family *fp;
  int k, found;

  for (fp = family_list; fp; fp = fp->next) {  /* ...loop over families */
    for (k = found = 0; k < nel; k++)          /* ...loop over elements */
      if (inFamily (U+k, fp))
	klist [found++] = k + 1;
    if (found) {                               /* ...show the members   */
      printf ("Family %3d: ", fp->id);
      for (k = 0; k < found; k++) 
	printf ("%d ", klist[k]);
      putchar ('\n');
    }
  }
  free (klist);
  return;
}
    
int showFamilies (void)
{
  Family *fp    = family_list;
  int     count = 0;

  if (use_families) {
    printf("Printing Family system...\n");
    for (fp = family_list; fp; fp = fp->next, count++) 
      printf("Family number %3d, members = %3d\n", count, fp->members);
    printf("Total number of families : %d\n", count);
  } else
    puts ("Family system is not active");

  return;
}

#endif

/* ------------------------------------------------------------------------ * 
 * This is the detection function that determines whether an element is in  *
 * a given family.  It checks the order of the element, then the element    *
 * perimeter to see if it matches.  To be in a family, an element must be   *
 * a simple shift of the family geometry.                                   * 
 * ------------------------------------------------------------------------ */

static int inFamily (Element *u, Family *f)
{
  int ok = 0;

  if (use_families) { 
    if (u->nr == f->nr && u->ns == f->ns) {
      int nb = (u->nr + u->ns - 2) << 1;
      double shift, dx, dy;
      const  double tol = dparam("TOLFAM");
      
      tempVector (fxb, nb); tempVector (fyb, nb);
      tempVector (uxb, nb); tempVector (uyb, nb);

      /* Get the (x,y)-coordinates of the perimeter */

      dgathr (nb, *f->parent->xmesh, f->parent->emap, fxb);
      dgathr (nb, *f->parent->ymesh, f->parent->emap, fyb);
      dgathr (nb, *   u     ->xmesh,    u     ->emap, uxb);
      dgathr (nb, *   u     ->ymesh,    u     ->emap, uyb);

      /* Shift the elements on top of each other */

      shift = *fxb - *uxb; dsadd (nb, shift, uxb, 1, uxb, 1);
      shift = *fyb - *uyb; dsadd (nb, shift, uyb, 1, uyb, 1);

      /* Check the difference between the two */

      dvsub (nb, fxb, 1, uxb, 1, uxb, 1);
      dvsub (nb, fyb, 1, uyb, 1, uyb, 1);

      dx = ddot(nb, uxb, 1, uxb, 1);
      dy = ddot(nb, uyb, 1, uyb, 1);
      if (sqrt(dx + dy) < tol) ok = 1;

      freeVector (fxb); freeVector (fyb);
      freeVector (uxb); freeVector (uyb);
    }
  } else                         /* ...not using families */
    { if ( u->id == f->id &&
           u->nr == f->nr && 
	   u->ns == f->ns ) ok = 1; }         

  return ok;
}

