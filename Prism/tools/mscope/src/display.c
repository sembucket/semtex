/*
 * DISPLAY routines for Mscope
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <stdarg.h>

#include "mscope.h"
#include "grammar.h"
#include "pl/pl.h"

extern Domain Geometry;

void DoDevice ()
{
  static char *help = 
    "x11             render to an X-window\n"
    "psfile <fname>  render to a PostScript file\n"
    "gif <fname>     render to a GIF file\n"
    "nil             render to the nil device (no graphics)\n";

  char *p = strtok(NULL,";\n");

  if (!p) {
    fprintf (stderr, "usage: dev <type> [options]\n");
    fprintf (stderr, "%s", help);
  } else {
    pl_device(p);
  }
}

/* Get limits and set them */

void GetLimits ()
{
  char *p;
  double limits [4];
  int n = 0;

  while (n < 4 && (p = strtok (NULL, " ")))
    limits[n++] = parse(p);
  
  if (n == 4)
    SetLimits (NULL, limits[0], limits[1], limits[2], limits[3]);
  else
    fputs ("usage: limits xmin xmax ymin ymax\n", stderr);

  return;
}

/* ------------------------------------------------------------------------- *
 * SetLimits                                                                 *
 *                                                                           *
 * These functions change the limits of the display.  If you call SetLimits  *
 * with gp != NULL, then it will calculate the limits for you based on the   *
 * problem geometry.  Otherwise, it uses the given values of (min,max) to    *
 * change the display limits.                                                *
 *                                                                           *
 * ZoomIn     Prompts for a box to zoom on.                                  *
 * ZoomOut    Resets the limits for the current geometry.                    *
 * ------------------------------------------------------------------------- */

void SetLimits (Domain *gp, ...)
{
  double dx, dy, scale;
  double xmin, xmax, ymin, ymax;

  if (gp) {
    dx    = dparam("XMAX") - (xmin = dparam("XMIN"));
    dy    = dparam("YMAX") - (ymin = dparam("YMIN"));
    scale = MAX(dx,dy);

    scale *= 1.05;   /* Leave some room around the edges */

    xmax  = xmin + 0.5 * ( dx + scale );
    xmin  = xmin + 0.5 * ( dx - scale );
    ymax  = ymin + 0.5 * ( dy + scale );
    ymin  = ymin + 0.5 * ( dy - scale );
  } else {
    va_list ap;
    va_start (ap, gp);
    xmin = (double) va_arg (ap, double);
    xmax = (double) va_arg (ap, double);
    ymin = (double) va_arg (ap, double);
    ymax = (double) va_arg (ap, double);
    va_end (ap);
  }

  pl_limits (xmin, xmax, ymin, ymax);

  dparam_set("mscope_xmin", xmin);
  dparam_set("mscope_xmax", xmax);
  dparam_set("mscope_ymin", ymin);
  dparam_set("mscope_ymax", ymax);
}

void ZoomIn ( void )
{
  int    k;
  double scale;
  double x1, x2, y1, y2;

  puts ("Enter points at opposite corners of the zoom box:");

  k = 0; while (k != 'e') pl_curs (&x1, &y1, &k);
  k = 0; while (k != 'e') pl_curs (&x2, &y2, &k);

  scale = MAX (fabs(x1-x2),fabs(y1-y2));
  x1    = (x1 + x2 - scale) * 0.5;
  x2    =  x1 + scale;
  y1    = (y1 + y2 - scale) * 0.5;
  y2    =  y1 + scale;

  pl_limits (x1, x2, y1, y2);
}

void DoErase ()
{
  pl_graphics();
  pl_erase   ();
  pl_gflush  ();
  pl_alpha   ();
}

void DoZoom() {
  if (Domain_require())
    ZoomIn();
}

void DoUnzoom() {
  if (Domain_require())
    SetLimits(&Geometry);
}

void DoRelocate()
{
  char *p = strtok(NULL, ";\n");

  if (!p) {
    fputs ("usage: relocate x0 y0\n", stderr);
  } else {
    char xexp[32];
    char yexp[32];

    if (sscanf(p,"%s%s", xexp, yexp) != 2)
      fprintf (stderr, "relocate: bad expression -- %s\n", p);

    pl_relocate (scalar(xexp), scalar(yexp));
  }
}
    
void DoDraw()
{
  char *p = strtok(NULL, ";\n");

  if (!p) {
    fputs ("usage: draw x0 y0\n", stderr);
  } else {
    char xexp[32];
    char yexp[32];

    if (sscanf(p,"%s%s", xexp, yexp) != 2)
      fprintf (stderr, "draw: bad expression -- %s\n", p);

    pl_draw  (scalar(xexp), scalar(yexp));
    pl_gflush();
  }
}

void DoLtype()
{
  char *p = strtok(NULL, ";\n");
  static char *usage = 
    "usage: ltype <type>, where <type> = \n"
    " 0 solid\n"
    " 1 dotted\n"
    " 2 dashed\n";

  if (!p)
    fputs (usage, stderr);
  else
    pl_ltype (atoi(p));
}

void DoLabel()
{
  char *p = strtok(NULL, ";\n");

  if (!p) {
    fputs ("usage: label <string>\n", stderr);
  } else {
    pl_label (p);
    pl_gflush();
  }
}
