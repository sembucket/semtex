/*
 * A simple xy plotting library
 *
 * Most of the function in PL are implemented in the driver.  This file only
 * contains the initialization code, and high-level drawing functions like
 * vector plots and contour plots.
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "pl.h"
#include "config.h"

double pl_xmin, pl_xmax;
double pl_ymin, pl_ymax;

/* ------------------------------------------------------------------------- */

void pl_init () {

  pl_xmin = pl_ymin = 0.;
  pl_xmax = pl_ymax = 1.;

  pl_driver_init();
}

void pl_exit () {
  pl_driver_exit();
}

/* ------------------------------------------------------------------------- */

static double *levs;     /* levels for contour plotting */
static int    nlevs;

void pl_levels (int n, double *x)
{
  if (levs) free(levs);

  nlevs = n;
  levs  = (double*) malloc(n*sizeof(double));
  memcpy(levs, x, n*sizeof(double));
}

/*
 * The values of icase correspond to the following lines:
 * icase == 1  => no line
 * icase == 7  => 2 and 8, or 5 and 3
 *
 *		       v10	         v11
 *			|---5----6----8---|
 *			|  5     6     8  |
 *			| 5      6      8 |
 *			|5       6       8|
 *			4444444444444444444
 *			|2       6       3|
 *			| 2      6      3 |
 *			|  2     6     3  |
 *			|___2____6____3___|
 *		       v00		 v01
 *		      (x,y)
 */

void pl_contour (int nx, int ny, double *xp, double *yp, double *up)
{
  double max_val;                 /* maximum allowed value */
  double no_val;                  /* value signifying missing data */
  double *vptr0, *vptr1;          /* pointers to rows in data */
  double *lptr;                   /* pointer to levels */
  double val;                     /* value of current level */
  double valmin, valmax;          /* attempt to limit level search */
  double v00, v01, v10, v11;      /* values in a cell */
  double x00, x01, x10, x11;      /* coordinates in a cell */
  double y00, y01, y10, y11;
  double xx0, xx1;                /* interpolated position */
  double yy0, yy1;
  
  int icase;                      /* which case to deal with */
  int i, j;                       /* counters in i and j */
  
  if (nlevs==0)
    fprintf (stderr, "pl: no levels for contour plotting!\n");
  else {

    max_val = FLT_MAX;
    no_val  = FLT_MAX;

    for (j = 0; j < ny-1; j++) {
      vptr0 = & up[j*nx];
      vptr1 = & up[(j+1)*nx];
      v01   = *vptr0++;
      v11   = *vptr1++;
      for (i = 0; i < nx-1; i++) {

	x00   = xp[j*nx + i];
	x01   = xp[j*nx + i + 1];
	x10   = xp[(j+1)*nx + i];
	x11   = xp[(j+1)*nx + i + 1];

	y00   = yp[j*ny + i];
	y01   = yp[j*ny + i + 1];
	y10   = yp[(j+1)*ny + i];
	y11   = yp[(j+1)*ny + i + 1];

	valmax = valmin = v00 = v01;

	if ((v01 = *vptr0++) < valmin) {
	  valmin = v01;
	} else if (v01 > valmax) {
	  valmax = v01;
	}
	if ((v10 = v11) < valmin) {
	  valmin = v10;
	} else if (v10 > valmax) {
	  valmax = v10;
	}
	if ((v11 = *vptr1++) != v11) { /* its' NaN */
	  fprintf(stderr,"pl: data value of %g at (%g,%g) is too large\n",
		  v11, x11, y11);
	  v11 = max_val;
	}

	if (v00 == no_val || v01 == no_val || v10 == no_val || v11 == no_val)
	  continue;

	if (v11 < valmin) {
	  valmin = v11;
	} else if (v11 > valmax) {
	  valmax = v11;
	}

	for (lptr = levs;
	     lptr < levs + nlevs && *lptr < valmin; lptr++) continue;

	while (lptr < levs + nlevs && *lptr <= valmax) {
	  val   = *lptr++;
	  icase = 1;

	  if (val > v00) icase++;
	  if (val > v01) icase += 2;
	  if (val > v10) icase += 4;
	  if (val > v11) icase  = 9 - icase;
	  
	  switch (icase) {
	  case 1:
	    continue;             /* continue to next level */
	  case 2:
	    xx0 = x00 + (x01-x00)*(val-v00)/(v01-v00);
	    yy0 = y00 + (y01-y00)*(val-v00)/(v01-v00);
	    xx1 = x00 + (x10-x00)*(val-v00)/(v10-v00);
	    yy1 = y00 + (y10-y00)*(val-v00)/(v10-v00);

	    break;
	  case 3:
	    xx0 = x00 + (x01-x00)*(val-v00)/(v01-v00);
	    yy0 = y00 + (y01-y00)*(val-v00)/(v01-v00);
	    xx1 = x01 + (x11-x01)*(val-v01)/(v11-v01);
	    yy1 = y01 + (y11-y01)*(val-v01)/(v11-v01);

	    break;
	  case 4:
	    xx0 = x00 + (x10-x00)*(val-v00)/(v10-v00);
	    yy0 = y00 + (y10-y00)*(val-v00)/(v10-v00);
	    xx1 = x01 + (x11-x01)*(val-v01)/(v11-v01);
	    yy1 = y01 + (y11-y01)*(val-v01)/(v11-v01);

	    break;
	  case 5:
	    xx0 = x00 + (x10-x00)*(val-v00)/(v10-v00);
	    yy0 = y00 + (y10-y00)*(val-v00)/(v10-v00);
	    xx1 = x10 + (x11-x10)*(val-v10)/(v11-v10);
	    yy1 = y10 + (y11-y10)*(val-v10)/(v11-v10);
	    
	    break;
	  case 6:
	    xx0 = x00 + (x01-x00)*(val-v00)/(v01-v00);
	    yy0 = y00 + (y01-y00)*(val-v00)/(v01-v00);
	    xx1 = x10 + (x11-x10)*(val-v10)/(v11-v10);
	    yy1 = y10 + (y11-y10)*(val-v10)/(v11-v10);

	    break;
	  case 7:			/* a saddle */

	    if(v00 > v01) {		/* draw cases 3 and 5 */

	      xx0 = x00 + (x01-x00)*(val-v00)/(v01-v00);
	      yy0 = y00 + (y01-y00)*(val-v00)/(v01-v00);
	      xx1 = x01 + (x11-x01)*(val-v01)/(v11-v01);
	      yy1 = y01 + (y11-y01)*(val-v01)/(v11-v01);
	      pl_line(xx0,yy0,xx1,yy1);
	      xx0 = x00 + (x10-x00)*(val-v00)/(v10-v00);
	      yy0 = y00 + (y10-y00)*(val-v00)/(v10-v00);
	      xx1 = x10 + (x11-x10)*(val-v10)/(v11-v10);
	      yy1 = y10 + (y11-y10)*(val-v10)/(v11-v10);

	    } else {			/* draw cases 2 and 8 */

	      xx0 = x00 + (x01-x00)*(val-v00)/(v01-v00);
	      yy0 = y00 + (y01-y00)*(val-v00)/(v01-v00);
	      xx1 = x00 + (x10-x00)*(val-v00)/(v10-v00);
	      yy1 = y00 + (y10-y00)*(val-v00)/(v10-v00);
	      pl_line(xx0,yy0,xx1,yy1);
	      xx0 = x10 + (x11-x10)*(val-v10)/(v11-v10);
	      yy0 = y10 + (y11-y10)*(val-v10)/(v11-v10);
	      xx1 = x01 + (x11-x01)*(val-v01)/(v11-v01);
	      yy1 = y01 + (y11-y01)*(val-v01)/(v11-v01);
	    }
	    break;
	  case 8:
	    xx0 = x10 + (x11-x10)*(val-v10)/(v11-v10);
	    yy0 = y10 + (y11-y10)*(val-v10)/(v11-v10);
	    xx1 = x01 + (x11-x01)*(val-v01)/(v11-v01);
	    yy1 = y01 + (y11-y01)*(val-v01)/(v11-v01);

	    break;
	  default:
	    fprintf(stderr,"pl: illegal value of icase %d\n",icase);
	    break;
	  }
	  pl_line(xx0,yy0,xx1,yy1);
	}
      }
    }
  }
}

/* ------------------------------------------------------------------------- */

static void vector (double x0, double y0, double len, double angle) 
{
  double c   = cos(angle);
  double s   = sin(angle);

  double xx0 = x0 - 0.5*len*c;  /*                 (xx2,yy2) \              */
  double yy0 = y0 - 0.5*len*s;  /*                            \             */
  double xx1 = x0 + 0.5*len*c;  /*                             \            */
  double yy1 = y0 + 0.5*len*s;  /*         ----------***---------           */
  double xx2 = x0 - 0.2*len*s;  /*  (xx0,yy0)      (x0,y0)     /  (xx1,yy1) */
  double yy2 = y0 + 0.2*len*c;  /*                            /             */
  double xx3 = x0 + 0.2*len*s;  /*                           /              */
  double yy3 = y0 - 0.2*len*c;  /*                 (xx3,yy3)                */

  pl_relocate (xx0, yy0);
  pl_draw     (xx1, yy1);
  pl_relocate (xx2, yy2);
  pl_draw     (xx1, yy1);
  pl_draw     (xx3, yy3);
}

void pl_vectors (int nx, int ny, double *x, double *y, double scale,
		 double *u, double *v)
{
  int i, j;
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      double xp    = x[j*nx + i];
      double yp    = y[j*nx + i];
      double up    = u[j*nx + i] + FLT_MIN;
      double vp    = v[j*nx + i] + FLT_MIN;
      double len   = hypot(up, vp);
      double angle = atan2(vp, up);
      vector (xp, yp, len*scale, angle);
    }
  }
}

/* ------------------------------------------------------------------------- *
 * Shade the area inside a curve.  This mimics a solid fill with closely     *
 * spaced lines in the current drawing color.  It looks nice for technical   *
 * drawings.                                                                 *
 * ------------------------------------------------------------------------- */

#define NCROSS 500

/* Rotate a set of positions by -alpha degrees */

static void rotate (int n, double alpha, double *a, double *b)
{
  const double cosa = cos(alpha);
  const double sina = sin(alpha);
  int i;

  for (i = 0; i < n; i++) {
    double tmp = cosa*a[i] + sina*b[i];
    b[i] = cosa*b[i] - sina*a[i];
    a[i] = tmp;
  }
}

/* Draw a line between (x0,y) and (x1,y) in rotated coordinates */

static void draw_seg (double x0, double x1, double y, double alpha) 
{
  const double cosa = cos(alpha);
  const double sina = sin(alpha);

  pl_relocate (x0*cosa - y*sina, y*cosa + x0*sina);
  pl_draw     (x1*cosa - y*sina, y*cosa + x1*sina);
}
 

#define TINY 1.e-5
#define LN2I 1.442695022		/* 1/ln(e) */
 
static void sort_dbl (int dimen, double *vec)
{
  double temp;
  int i,j,m,n;				/* counters */
  int lognb2;				/* (int)(log_2(dimen)) */

  if (dimen <= 0) return;

  lognb2 = log((double)dimen)*LN2I + TINY;	/* ~ log_2(dimen) */
  m = dimen;
  for (n = 0; n < lognb2; n++) {
    m /= 2;
    for (j = m; j < dimen; j++) {
      i = j - m;
      temp = vec[j];
      while (i >= 0 && temp < vec[i]) {
	vec[i+m] = vec[i];
	i -= m;
      }
      vec[i+m] = temp;
    }
  }
}
 
#ifdef DRIVER_sm

void pl_shade (int n, double *x, double *y, int density)
{
  int i;

  float *fx = (float*) malloc(n*sizeof(float));
  float *fy = (float*) malloc(n*sizeof(float));

  for (i = 0; i < n; i++) {
    fx[i] = (float) x[i];
    fy[i] = (float) y[i];
  }

  sm_shade (density, fx, fy, n);

  free(fx);
  free(fy);
}

#else

void pl_shade (int n, double *x, double *y, int delta) 
{
  const double angle = 45*(M_PI/180.);

  double dy;   /* line spacing in y units */
  double yval; /* current (rotated) y value */
  int i, j, i1;

  double crossing[NCROSS];     /* where line segments cross curve */
  double cx[4];                /* corners of box */
  double cy[4];

  cx[0] = pl_xmin;    cy[0] = pl_ymin;
  cx[1] = pl_xmax;    cy[1] = pl_ymin;
  cx[2] = pl_xmin;    cy[2] = pl_ymax;
  cx[3] = pl_xmax;    cy[3] = pl_ymax;

  rotate(n, angle,  x,  y);
  rotate(4, angle, cx, cy);

  sort_dbl (4, cx);
  sort_dbl (4, cy);

  dy = hypot(pl_xmax-pl_xmin,pl_ymax-pl_ymin)/delta;

  /* Now draw the lines */

  for (yval = cy[0]; yval < cy[3]; yval += dy) {
    for (i = j = 0; i < n; i++) {
      i1 = (i == n-1) ? 0 : i + 1;
      if (i1 != 0 && (y[i1]-yval)==0.0) {
	continue;  /* dont' find it twice */
      }
      if ((y[i] - yval)*(y[i1]-yval) <= 0.0) { /* a crossing */
	if ((yval-y[i]) == 0.0) {
	  crossing[j++] = x[i];    /* avoid division by zero */
	} else {
	  crossing[j++] = x[i] + 
	    (x[i1]-x[i])/(y[i1] - y[i])*(yval-y[i]);
	}
	if (j > NCROSS) {
	  fprintf (stderr, "pl: too many crossings!\n");
	  break;
	}
      }
    }
    if (j > 0) {                  /* lines to draw */
      sort_dbl(j, crossing);
      for (i = 0; i < j-1; i += 2)
	draw_seg(crossing[i], crossing[i+1], yval, angle);
      if (j%2 == 1) {             /* an incomplete line */
	draw_seg(crossing[j-1],cx[3],yval, angle);
      }
    }
  }

  rotate(n, -angle, x, y);
}
#endif

