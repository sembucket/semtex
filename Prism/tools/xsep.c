/*
 * Compute separation points
 *
 * In general, "separation" occurs in a flow when the flow ceases to make
 * contact with a solid boundary.  Mathematically, this is indicated when
 * the normal shear stress goes to zero, i.e. tau_n(xs,ys) = 0.
 *
 * This program finds separation points by searching along a line (y=const)
 * to find the value of xs such that tau_n(xs,y) = 0.  It's important to have
 * a good initial guess for these points! Newton's method is used to refine 
 * the initial guess and solve for x.
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"

char *prog   = "xsep";
char *author = "Ron Henderson";
char *rcsid  = "$Id$";
char *usage  = "usage: xsep session.rea session.fld\n";
char *help   =
"This program computes separation points, tau_n(xs,ys)=0, given an initial  \n"
"guess for where they occur.  During the search, y is held fixed and x is   \n"
"varied by applying Newton's method.\n";

/* ------------------------------------------------------------------------- */

main (int argc, char *argv[])
{
  FILE *fp;
  FieldFile *f;

  Field *U;
  Field *V;
  Field *tau;
  Field *work;

  Bedge *Ubc, *Walls, *bc;

  if (argc != 3) {
    fputs (usage, stderr);
    exit (-1);
  }

  /* ---------- Read the input file ---------- */

  fp = fopen(argv[1], "r");
  
  speclib_init();
  ReadParams(fp);
  U    = ReadMesh(fp);
  V    = Field_dup(U); FIELD_TYPE(V) = 'v';
  tau  = Field_dup(U);
  work = Field_dup(U);
  Ubc  = ReadBCs (fp, 0, U);
  fclose (fp);

  /* ---------- Read the field file ---------- */

  fp = fopen(argv[2], "r");
  f  = (FieldFile*) calloc(1, sizeof(FieldFile));
  FieldFile_read (f, fp);
  FieldFile_get  (f, U);
  FieldFile_get  (f, V);
  FieldFile_free (f);
  fclose (fp);

  /* Compute tau_{1,2}, the only component of the stress tensor that     *
   * contributes to the fluid drag along the walls.  Just for reference, *
   *                                                                     *
   *    tau_{1,2} = (d u_1/d x_2 + d u_2/d x_1) = ( dU/dy + dV/dx )      *
   *                                                                     *
   * Now locate the separation points by finding tau(x,y) = 0.  For the  *
   * Newton iteration we need both tau and \partial_x tau, which is      *
   * computed and stored in the Field "work".                            */

  Field_dy   (U, tau);
  Field_dx   (V, work);
  Field_axpy (1., work, tau);
  Field_dx   (tau, work);

  printf ("%s: type Control-C to kill me\n\n", prog);
  while (1) {
    double x, y;
    printf ("Enter initial guess [x,y]: ");
    scanf  ("%lf%lf", &x, &y);
    locate_xsep (tau, work, x, y);
  }

  return 0;
}

/* Locate a separation point by using Newton's method to solve for  *
 * the point x where tau(x,y) = 0.                                  */

int locate_xsep (Field *tau, Field *tau_x, double xo, double y)
{
  const double tol   = dparam("TOLABS");
  const int max_iter = 10;

  Probe *probe;
  double x, dx;
  double f, df;
  int iter = 0;

  printf ("Newton iteration, tolerance = %g\n", tol);

  x = xo;
  probe = Probe_alloc (tau, PROBE_XP, x, y);
  do {
    Probe_move(probe, x, y);

    f  = Probe_eval (probe, tau);
    df = Probe_eval(probe, tau_x);
    dx = f/df;

    printf ("iter %d: x = %g, dx = %g, f = %g, df = %g\n", iter, x, dx, f, df);

    x -= dx;
  } 
  while (fabs(dx) > tol && ++iter < max_iter);

  Probe_free (probe);
  return 0;
}
