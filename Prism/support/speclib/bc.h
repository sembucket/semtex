#ifndef BC_H
#define BC_H

/* BCs
 *
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 * 
 * ------------------------------------------------------------------------- */

#include "speclib/edge.h"
#include "speclib/element.h"
#include "speclib/field.h"

typedef enum {
  NEUMANN,
  DIRICHLET,
  UNKNOWN
} BC_type;

typedef struct bedge {              /* .... BOUNDARY EDGE Definition ..... */
  int             id         ;      /* ID number                           */
  char            type       ;      /* Type (Flux, Temperature, etc.)      */
  union {                           /* Boundary condition info:            */
    double        *value     ;      /*    - value                          */
    char          *function  ;      /*    - function                       */
  } bc;                             /*                                     */
  struct element  *elmt      ;      /* Corresponding element               */
  struct edge     *edge      ;      /* Corresponding edge                  */
  struct bedge    *next      ;      /* Pointer to the next one             */
} Bedge;

/* Prototypes */

int      BC_init          (void);
Bedge*   BC_alloc         (void);
void     BC_free          (Bedge *list);
int      BC_count         (Bedge *list);
Bedge*   BC_get           (Bedge *list, char type);
Bedge*   BC_dup           (Bedge *list);
void     BC_learn         (Bedge *list, Field *u);
Bedge*   BC_make          (char type, Element *elmt, Edge *edge, ...);

/* The library maintains a map between char tags and BC_type's.  You can add *
 * user-defined tags, override default tags, or check a tag's type using the *
 * functions below.   For example, to define tag 'A' as a DIRICHLET boundary *
 * condition you would use BC_defType():                                     *
 *                                                                           *
 *    BC_defType ('A', DIRICHLET);                                           */

BC_type  BC_getType       (char c);
BC_type  BC_defType       (char c, BC_type type);

#endif













