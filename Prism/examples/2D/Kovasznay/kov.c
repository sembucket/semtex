/*
 * Kovasznay Flow - an exact solution to the Navier-Stokes Equations  
 * 
 * $Id$
 *
 * Author: R. D. Henderson
 */

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "prism/prism.h"

#ifndef   KOVXY
#  define KOVXY
#endif

#ifdef    KOVXZ
#  undef  KOVXZ
#endif

static struct {
  struct {
    char  type;
    char* solution;
  } u[2];
  struct {
    char  type;
    char* solution;
  } p;
} Kovasznay = {
  'u',    "1-exp(lambda*x)*cos(2*PI*y)",
  'v',    "lambda*exp(lambda*x)*sin(2*PI*y)/(2*PI)",
  'p',    "(lambda^2/(2*PI)-2*PI)*exp(lambda*x)*sin(2*PI*y)"
};

static Field *u_exact[DIM];
static Field *p_exact;

/* ----------------------------------------------------------------------- */

void Kovasznay_analyze (Domain *domain)
{
  const int    step = iparam("STEP");
  const double time = dparam("TIME");

  Field *u[3];
  Field *p;
  Field *error;

  int i;

  u[0] = domain->U;
  u[1] = domain->V;
  u[2] = domain->W;
  p    = domain->P;

  /* initialize the exact solution */

  if (u_exact[0] == (Field*) NULL) {
    for (i = 0; i < DIM; i++) {
      u_exact[i] = Field_dup(u[i]);
      FIELD_TYPE(u_exact[i]) = Kovasznay.u[i].type;
      Field_set (u_exact[i], Kovasznay.u[i].solution);
    }

    p_exact = Field_dup(p);
    FIELD_TYPE(p_exact) = Kovasznay.p.type;
    Field_set (p_exact, Kovasznay.p.solution);
  }

  error = Field_dup(u[0]);
  FIELD_TYPE(error) = 'e';

  /* compute velocity error */

  Field_copy (u_exact[0], error);
  Field_axpy (-1., u[0],  error);
  printf ("%#14.7g %#14.7g", Field_L2(error), Field_H1(error));

  Field_copy (u_exact[1], error);
  Field_axpy (-1., u[1],  error);
  printf ("%#14.7g %#14.7g", Field_L2(error), Field_H1(error));
  putchar('\n');

  Field_free(error);
}

/* ------------------------------------------------------------------------- */

int main (int argc, char *argv[]) 
{
  Domain *domain;
  double  t0, tn;

  Prism_init(&argc,&argv);
  parse_args( argc, argv);

  domain = Domain_alloc(argv[argc-1]);
  t0     = dparam("TIME_0");
  tn     = dparam("TIME_N");

  Navier_Stokes (domain,t0,tn);
  PostProcess   (domain);

  printf ("u min,max = %g %g\n", Field_min(domain->U), Field_max(domain->U));
  printf ("v min,max = %g %g\n", Field_min(domain->V), Field_max(domain->V));

  Kovasznay_analyze(domain);

  Domain_free(domain);
  Prism_exit();
}
