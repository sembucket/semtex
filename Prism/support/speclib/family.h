#ifndef FAMILY_H
#define FAMILY_H

/* Family
 *
 * $Id$
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * ----------------------------------------------------------------------- */

#include "speclib/element.h"
#include "speclib/field.h"

typedef struct family {             /* ........ FAMILY Definition ........ */
  int            id             ;   /* ID number of the parent element     */
  int            members        ;   /* Number of members                   */
  int            set            ;   /* Set flag (on = 1, off = 0)          */
  int            nr, ns         ;   /*                                     */
  double         **xr, **xs     ;   /* ----------------------------------- */
  double         **yr, **ys     ;   /*                                     */
  double         **jac          ;   /*       F A M I L Y   D A T A         */
  double         **rx, **ry     ;   /*             (Geometry)              */
  double         **sx, **sy     ;   /*                                     */
  double         **mass         ;   /* ----------------------------------- */
  struct element *parent        ;   /* Defining element for the family     */
  struct family  *next          ;   /* Pointer to the next family          */
} Family;


/* Prototypes */

Family*  Family_get       (const Element *elmt);
Family*  Family_create    (Element *elmt);
int      Family_count     (const Field *u);
int      Family_set       (Family  *f);
int      Family_add       (Family  *f, Element *elmt);
void     Family_reset     (void);
void     Family_disable   (void);
void     Family_destroy   (void);

#endif
