/* 
 * Family system
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
 * Here is a list of the methods for working with element Families:
 *
 *     get             - Get the family for a particular element
 *     create          - Create a new family
 *     add             - Add an element to the given family         
 *
 *     set             - Mark a family as having been processed
 *     reset           - Reset all families                   
 *     count           - Count the number of element families in a mesh
 *
 *     destroy         - Free all family information
 *                                                                          
 * The functions set() and reset() allow you to perform a "family loop."  
 * The lists and effects of families are otherwise invisible outside of this 
 * file (This is flaw --- the use of a static family list causes problems
 * in some applications).
 *                                                                          
 * A call to the function Family_disable() turns off the use of families.
 *
 * A call to Family_destroy() frees all Family information and resets the
 * system.  This should be used before loading a new mesh, for example.
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
#include "veclib/veclib.h"
#include "speclib/speclib.h"

#define GEOMEMBERS  10    /* # of geometry factors in a family */

/* Externals */

static int use_families = 1;
static Family *family_list;

/* ------------------------------------------------------------------------ */

/* The following function should be called once before any  *
 * pre-processing to disable the use of element families.   */

void Family_disable (void) { use_families = 0; return; }

/* ------------------------------------------------------------------------ */

void Family_destroy (void)
{
  Family *fam = family_list;

  while (fam) {
    Family *next = fam->next;
    free(fam);
    fam = next;
  }

  family_list = NULL;
}

/* ------------------------------------------------------------------------ */

static Family *append(Family *list, Family *f)
{
  Family *t;
  Family *p = 0;
  int i = 0;
  
  for (t = list; t; t = t->next) 
    { p = t; ++i; }

  if (p) {
    p->next = f;
    f->id   = i;
  } else
    list = f;

  return list;
}

Family *Family_create (Element *elmt)
{
  Family *fam = (Family*) calloc (1, sizeof(Family));
  
  /* Set family characteristics */
  
  fam->nr       = elmt->nr;
  fam->ns       = elmt->ns;
  fam->members  = 1; 
  fam->set      = 0;

  /* Save the geometry data for this family */

  fam->xr   = elmt->xr;
  fam->xs   = elmt->xs;
  fam->yr   = elmt->yr;
  fam->ys   = elmt->ys;
  fam->jac  = elmt->jac;
  fam->rx   = elmt->rx;
  fam->ry   = elmt->ry;
  fam->sx   = elmt->sx;
  fam->sy   = elmt->sy;
  fam->mass = elmt->mass;

  /* Add to the family list */

  fam->parent  = elmt;
  family_list  = append(family_list, fam);

  if (!use_families) fam->id = elmt->id;

  return fam;
}

/* This just adds the given element to an existing family.  Returns *
 * the new number of members in the family.                         */

int Family_add (Family *f, Element *elmt)
{
  elmt->xr   = f->xr;
  elmt->xs   = f->xs;
  elmt->yr   = f->yr;
  elmt->ys   = f->ys;
  elmt->jac  = f->jac;
  elmt->rx   = f->rx;
  elmt->ry   = f->ry;
  elmt->sx   = f->sx;
  elmt->sy   = f->sy;
  elmt->mass = f->mass;

  return ++(f->members);
}

/* ------------------------------------------------------------------------ * 
 * This is the detection function that determines whether an element is in  *
 * a given family.  It checks the order of the element, then the element    *
 * perimeter to see if it matches.  To be in a family, an element must be   *
 * a simple shift of the family geometry.                                   * 
 * ------------------------------------------------------------------------ */

static int isInFamily (const Family *f, const Element *elmt)
{
  int ok = 0;

  if (use_families) { 
    if (elmt->nr == f->nr && elmt->ns == f->ns) {
      int nb = (elmt->nr + elmt->ns - 2) << 1;
      double shift, dx, dy;
      const double tol = dparam("TOLFAM");
      
      tempVector (fxb, nb); tempVector (fyb, nb);
      tempVector (uxb, nb); tempVector (uyb, nb);

      /* Get the (x,y)-coordinates of the perimeter */

      dgathr (nb, *f->parent->xmesh, f->parent->emap, fxb);
      dgathr (nb, *f->parent->ymesh, f->parent->emap, fyb);
      dgathr (nb, *elmt     ->xmesh, elmt     ->emap, uxb);
      dgathr (nb, *elmt     ->ymesh, elmt     ->emap, uyb);

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
  } else {                       /* ...not using families */
    if (elmt->id == f->id &&
	elmt->nr == f->nr && 
	elmt->ns == f->ns) 
      ok = 1; 
  }         

  return ok;
}

/* Returns a pointer to the family for a given element.  If the family    *
 * system is disabled, this function just matches ID numbers.  Otherwise  *
 * it matches the perimeter of the element and the family.                *
 *                                                                        *
 * Returns NULL if no family is found.                                    */

Family *Family_get (const Element *elmt)
{
  Family *f;
  for (f = family_list; f; f = f->next)
    if (isInFamily (f, elmt)) break;
  return f;
}

/* Mark this family as set...for doing loops over families. *
 * Returns the number of times the family has been set.     */

int Family_set (Family *fam) { return ++(fam->set); }

/* Turn the "set" flag off for all families */

void Family_reset (void)
{ 
  Family *f = family_list;
  while (f) 
    { f->set = 0; f = f->next; }
}

/* This returns the number of families in a linked element list. *
 * Can be less than the number of members of the family list.    */

int Family_count (const Field *U)
{
  const int nel = Field_count(U);
  int count = 0;
  int k;

  if (!use_families) return nel;

  Family_reset();
  for (k = 0; k < nel; k++) {
    Family *f = Family_get(&U[k]);
    if (f && f->set != 1) {
      count++;
      Family_set (f);
    }
  }

  return count;
}

