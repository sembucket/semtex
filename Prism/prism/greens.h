#ifndef GREENS_H
#define GREENS_H

/* $Id$ */
 
#include "prism/constants.h"
#include "speclib/field.h"

typedef struct gf_ {              /* ....... Green's Function ........ */
  int       order               ; /* Time-order of the current field   */
  Bedge    *Gbc                 ; /* Special boundary condition array  */
  Field    *basis               ; /* Basis velocity (U or W)           */
  Field    *Gv[DIM][_MAX_TORDER]; /* Green's function velocities       */
  Field    *Gp     [_MAX_TORDER]; /* Green's function pressure         */
  double    Fg     [_MAX_TORDER]; /* Green's function "force"          */
} GreensF;

#endif
