/*
 * PL driver based on GNU libplot
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <plot.h>
#include "pl.h"

/* ------------------------------------------------------------------------- */

static int    handle;               /* current plotter */
static int    devnum;               /* current device number */
static int    marker_type = M_DOT;  /* current marker type */
static double marker_size = 1.;     /* current marker size */
static double expand = 1.;          /* current expand ratio */
static double angle  = 0.;          /* current angle (for text?) */

extern double pl_xmin;              /* current user limits */
extern double pl_xmax;
extern double pl_ymin;
extern double pl_ymax;

#ifndef MIN
#define MIN(a,b)  ((a)<(b)?(a):(b))
#endif

/* ------------------------------------------------------------------------- */

static int nodev_handle = -1;

static int nodev_open ()  
{ 
  if (nodev_handle == -1) {
    FILE *null = fopen("/dev/null","w");
    nodev_handle = newpl("meta", stdin, null, stderr);
    selectpl(nodev_handle);
    openpl();
  }

  selectpl(handle = nodev_handle);
  return 0; 
}

static int nodev_close() { return 0; }

/* ------------------------------------------------------------------------- */

static int x11_handle = -1;

static int x11_open()  
{ 
  if (x11_handle == -1) {
    x11_handle = newpl("X", stdin, stdout, stderr);
    selectpl(x11_handle);
    openpl();
  }
   
  selectpl(handle = x11_handle);
  pl_limits(pl_xmin,pl_xmax,pl_ymin,pl_ymax);
  return 0; 
}

static int x11_close() { 
  flushpl();
  return 0; 
}

/* ------------------------------------------------------------------------- */

static int ps_open (char *args) 
{
  char fname[FILENAME_MAX];
  sscanf(args, "%*s%s", fname);

  handle = newpl("ps", stdin, fopen(fname,"w"), stderr);
  selectpl(handle);
  openpl();

  pl_limits(pl_xmin,pl_xmax,pl_ymin,pl_ymax);

  return 0;
}

static int ps_close() 
{ 
  int ps_handle = handle;

  flushpl();
  closepl();
  selectpl(handle=0);
  deletepl(ps_handle);

  return 0; 
}

/* ------------------------------------------------------------------------- */

static char gif_file[FILENAME_MAX];
static char gif_temp[FILENAME_MAX];

static int gif_open (char *args)  
{ 
  sscanf (args, "%*s%s", gif_file);
  sprintf(gif_temp, "ps_%d", rand());

  handle = newpl("ps", stdin, fopen(gif_temp,"w"), stderr);
  selectpl(handle);
  openpl();
  pl_limits(pl_xmin,pl_xmax,pl_ymin,pl_ymax);

  return 0; 
}

static int gif_close() 
{ 
  int gif_handle = handle;
  char buf[BUFSIZ];

  flushpl();
  closepl();
  selectpl(handle=0);
  deletepl(gif_handle);

  sprintf(buf,
	  "ghostscript -q -sOutputFile=- -sDEVICE=ppm - < %s "
	  "| ppmtogif > %s", gif_temp, gif_file);
  system(buf);
  unlink(gif_temp);

  return 0; 
}

/* ------------------------------------------------------------------------- *
 * Supported Devices                                                         *
 *                                                                           *
 * For a device to be accessed through PL, it must have an entry in the      *
 * following table.  All that's necessary is the name of the device, and     *
 * pointers to the functions that open and close it.  Some devices, like     *
 * the x11 device, cannot be closed, so that open() simply reactivates the   *
 * device.   Note that this automatically eliminates the persistence of      *
 * libplot X plotters, since they are only forked when closed.               *
 * ------------------------------------------------------------------------- */

typedef struct {
  char *name;          /* name of the device */
  int (*open)();       /* open the device */
  int (*close)();      /* close the device */
} device_t;

static device_t devTable[] = {
  {"nodevice",  nodev_open, nodev_close },
  {"x11",       x11_open,   x11_close   },
  {"postfile",  ps_open,    ps_close    },
  {"gif",       gif_open,   gif_close   },

  /* aliases */

  {"X",         x11_open,   x11_close   },
  {"null",      nodev_open, nodev_close },
  {"nil",       nodev_open, nodev_close },
  {"ps",        ps_open,    ps_close    },

  { 0, 0, 0 }
};

void pl_driver_init () {
  devnum = 0;
  devTable[devnum].open();
}

void pl_driver_exit () {
  pl_device("nodevice");
}

void pl_device (char *dev) 
{
  /* close the current device */

  if (devnum != 0) devTable[devnum].close();  

  /* Search the table for a name matching the requested device */

  for (devnum = 0; devTable[devnum].name != NULL; devnum++) {
    const char *target = devTable[devnum].name;
    if (strncmp(dev, target, strlen(target))==0) {
      devTable[devnum].open(dev);
      break;
    }
  }
}							     

/* ------------------------------------------------------------------------- */

void pl_erase() {
  erase();
  flushpl();
}

void pl_gflush() {
  flushpl();
}

void pl_label (const char *str) {
  alabel ('c', 'c', str);
  flushpl();
}

void pl_expand (double val)  {}

void pl_angle   (double val) {
  ftextangle(val);
}

void pl_limits (double xmin, double xmax, double ymin, double ymax) 
{
  double x0, y0, x1, y1;

  /* Save coordinates of the user-drawable window */

  pl_xmin = xmin;
  pl_xmax = xmax;
  pl_ymin = ymin;
  pl_ymax = ymax;

  /* Increase limits by a little bit to allow for axes, etc. */

  x0 = xmin - (xmax-xmin)*0.1;
  x1 = xmax + (xmax-xmin)*0.1;
  y0 = ymin - (ymax-ymin)*0.1;
  y1 = ymax + (ymax-ymin)*0.1;

  fspace(x0,y0,x1,y1);
  
  /* Reset default values for scaled quantities */

  marker_size = 0.02*MIN(xmax-xmin,ymax-ymin);
}

void pl_connect (int n, double *x, double *y) 
{
  fmove(*x,*y);
  while (--n) {
    fcont(*++x, *++y);
  }
  endpath();
  flushpl();
}

void pl_points  (int n, double *x, double *y) 
{
  while (n--)
    fmarker(*x++, *y++, marker_type, marker_size);
  flushpl();
}

void pl_dot () {
  fmarkerrel (0., 0., marker_type, marker_size);
  flushpl();
}

void pl_box()
{
  fmove(pl_xmin,pl_ymin);
  fcont(pl_xmax,pl_ymin);
  fcont(pl_xmax,pl_ymax);
  fcont(pl_xmin,pl_ymax);
  fcont(pl_xmin,pl_ymin);
  endpath();
  flushpl();
}

void pl_xlabel  (const char *str) {}
void pl_ylabel  (const char *str) {}

void pl_line (double x0, double y0, double x1, double y1) {
  fline(x0,y0,x1,y1);
  flushpl();
}

void pl_relocate(double x, double y) {
  fmove(x,y);
}

void pl_draw (double x, double y) {
  fcont(x,y);
}

void pl_ltype (int type) 
{
  switch (type) {
  case 0:
    linemod("solid");
    break;
  case 1:
    linemod("dotted");
    break;
  case 2:
    linemod("shortdashed");
    break;
  case 3:
    linemod("longdashed");
    break;
  case 4:
    linemod("dotdashed");
    break;
  default:
    linemod(NULL);
    break;
  }
}

void pl_ptype (int type) 
{
  switch(type) {
  case PL_DOT:
    marker_type = M_DOT;
    break;
  case PL_PLUS:
    marker_type = M_PLUS;
    break;
  case PL_ASTERISK:
    marker_type = M_ASTERISK;
    break;
  case PL_CIRCLE:
    marker_type = M_CIRCLE;
    break;
  case PL_CROSS:
    marker_type = M_CROSS;
    break;
  case PL_STAR:
    marker_type = M_STAR;
    break;
  case PL_SQUARE:
    marker_type = M_SQUARE;
    break;
  case PL_DIAMOND:
    marker_type = M_DIAMOND;
    break;
  case PL_TRIANGLE:
    marker_type = M_TRIANGLE;
    break;
  case PL_INVERTED_TRIANGLE:
    marker_type = M_INVERTED_TRIANGLE;
    break;
  case PL_FILLED_SQUARE:
    marker_type = M_FILLED_SQUARE;
    break;
  case PL_FILLED_DIAMOND:
    marker_type = M_FILLED_DIAMOND;
    break;
  case PL_FILLED_TRIANGLE:
    marker_type = M_FILLED_TRIANGLE;
    break;
  case PL_FILLED_INVERTED_TRIANGLE:
    marker_type = M_FILLED_INVERTED_TRIANGLE;
    break;
  case PL_FILLED_CIRCLE:
    marker_type = M_FILLED_CIRCLE;
    break;
  case PL_STARBURST:
    marker_type = M_STARBURST;
    break;
  case PL_FANCY_PLUS:
    marker_type = M_FANCY_PLUS;
    break;
  case PL_FANCY_CROSS:
    marker_type = M_FANCY_CROSS;
    break;
  case PL_FANCY_SQUARE:
    marker_type = M_FANCY_SQUARE;
    break;
  case PL_FANCY_DIAMOND:
    marker_type = M_FANCY_DIAMOND;
    break;
  case PL_FILLED_FANCY_SQUARE:
    marker_type = M_FILLED_FANCY_SQUARE;
    break;
  case PL_FILLED_FANCY_DIAMOND:
    marker_type = M_FILLED_FANCY_DIAMOND;
    break;
  case PL_HALF_FILLED_SQUARE:
    marker_type = M_HALF_FILLED_SQUARE;
    break;
  case PL_HALF_FILLED_DIAMOND:
    marker_type = M_HALF_FILLED_DIAMOND;
    break;
  case PL_HALF_FILLED_TRIANGLE:
    marker_type = M_HALF_FILLED_TRIANGLE;
    break;
  case PL_HALF_FILLED_INVERTED_TRIANGLE:
    marker_type = M_HALF_FILLED_INVERTED_TRIANGLE;
    break;
  case PL_HALF_FILLED_CIRCLE:
    marker_type = M_HALF_FILLED_CIRCLE;
    break;
  case PL_HALF_FILLED_FANCY_DIAMOND:
    marker_type = M_HALF_FILLED_FANCY_DIAMOND;
    break;
  case PL_OCTAGON:
    marker_type = M_OCTAGON;
    break;
  case PL_FILLED_OCTAGON:
    marker_type = M_FILLED_OCTAGON;
    break;
  default:
    marker_type = M_NONE;
    break;
  }
}

/* The following are obsolete */

void pl_curs    (double *x, double *y, int *key) {}
void pl_zoom    () {}
void pl_alpha()    {}
void pl_graphics() {}




