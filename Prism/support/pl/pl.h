#ifndef PL_H
#define PL_H

/* PL
 *
 * $Id$
 *
 * Author: R. D. Henderson
 *
 * PL is a simple 2D vector graphics library.
 * ------------------------------------------------------------------------- */

#include "pl/ptype.h"
#include "pl/ltype.h"
#include "pl/devices.h"

void pl_init     (void);
void pl_exit     (void);

void pl_device   (char *dev);
void pl_gflush   ();
void pl_alpha    ();
void pl_graphics ();

void pl_erase    ();
void pl_limits   (double xmin, double xmax, double ymin, double ymax);
void pl_angle    (double val);
void pl_expand   (double val);

void pl_ltype    (int type);
void pl_ptype    (int type);

void pl_relocate (double x,  double y);
void pl_draw     (double x,  double y);
void pl_line     (double x0, double x1, double y0, double y1);
void pl_dot      ();
void pl_points   (int n, double *x, double *y);
void pl_connect  (int n, double *x, double *y);
void pl_shade    (int n, double *x, double *y, int density);

void pl_label    (const char *str);
void pl_xlabel   (const char *str);
void pl_ylabel   (const char *str);
void pl_box      ();

void pl_levels   (int n, double *levels);
void pl_contour  (int nx, int ny, double *x, double *y, double *u);
void pl_vectors  (int nx, int ny, double *x, double *y, double scale,
		  double *u, double *v);

#endif
