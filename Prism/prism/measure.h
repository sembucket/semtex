#ifndef MEASURE_H
#define MEASURE_H

/* measure_t
 *
 * This data structure is used to take measurements of body forces and mode 
 * amplitudes in a 2D or 3D field.   A measure_t is allocated by giving a
 * pointer to a domain and a file name.  Each time you call the function
 * measure_analyze() it performs its calculations and writes the data to
 * the named output file.
 *
 * The output (measurements) file is simply a text file with a short header
 * describing the contents.
 *
 * Note that measure_t is simple a subset of Domain and not a stand-alone
 * data structure.
 *
 * Copyright (c) 1997 Ronald D. Henderson and Caltech
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "speclib/speclib.h"

/* macros to access measured quantities */

#define MEASURE_FP(m,i)     (m->force[i].pressure)
#define MEASURE_FV(m,i)     (m->force[i].viscous)
#define MEASURE_FORCE(m,i)  (MEASURE_FP(m,i) + MEASURE_FV(m,i))
#define MEASURE_ENERGY(m,k) (m->energy[(k)])
#define MEASURE_AREA(m)     (m->area)

typedef struct measure {
  FILE*  fp;         /* Output stream */
  Bedge* body;       /* List of boundary segments */
  double viscosity;  /* Fluid viscosity (for skin friction) */
  double area;       /* Domain size (for normalization) */

  struct {           /* Forces in each direction */
    double pressure; /*    pressure */
    double viscous;  /*    viscous */
  } force[3];
  double *energy;    /* Integrated kinetic energy */

  struct domain *d ; /* Computational domain */
} measure_t;

/* Prototypes */

measure_t *measure_alloc   (struct domain *d);
void       measure_free    (measure_t *m);
void       measure_analyze (measure_t *m);

#endif











