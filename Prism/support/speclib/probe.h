#ifndef PROBE_H
#define PROBE_H

/* Probe
 *
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 *
 * A Probe provides a method for evaluating a Field at an arbitrary point,
 * or at the discrete mesh point closest to a given position.
 * ------------------------------------------------------------------------- */

#include "speclib/field.h"

typedef enum {
  PROBE_XP,
  PROBE_INDEX
} location_t;

typedef struct probe {       /* ------------  Probe ------------- */
  double        x, y;        /* Coordinates in (x,y)-space        */
  location_t    type;        /* Position or Index                 */
  union {                    /* Location info                     */
    struct {
      double r, *hr;         /* Coordinates in (r,s)-space        */
      double s, *hs;         /*                                   */
    } mesh;
    struct {
      int    ir, is;         /* Coordinates in (i,j)-space        */
    } node;
  } location;
  Field*        geom;        /* Field that provides the geometry  */
  Element*      elmt;        /* Element the probe lies within     */
} Probe;


Probe *Probe_alloc (const Field *u, location_t type, double x, double y);
void   Probe_free  (Probe *p);

double Probe_eval (const Probe *p, const Field *u);

/* Reset the location of a probe.  This function returns 0 if everything  *
 * goes OK.  If the new point (x,y) cannot be located within the mesh,    *
 * the probe is left unchanged and the return value is non-zero.          */

int Probe_move (Probe *p, double x, double y);

#endif
