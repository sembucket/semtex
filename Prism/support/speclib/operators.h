#ifndef OPERATORS_H
#define OPERATORS_H

/* Operators
 *
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 *
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdarg.h>

#include "speclib/element.h"
#include "speclib/field.h"
#include "speclib/matrix.h"

void coef   (int np);
void getops (int np, double* *zp, double* *wp, double* **dp, double* **dpt);

void geofac      (Element *U, Matrix_IP *G);
void Helmholtz   (Field   *U, BSystem *B, double *in, double *out);
void HelmholtzSC (Element *U, BSystem *B, 
		  double **A_bb, double **A_bi, double **A_ii);

#endif
