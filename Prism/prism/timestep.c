/*
 * Time Stepping and Integration Routines
 *
 * $Id$
 * --------------------------------------------------------------------- */

#include <math.h>

#include "prism/prism.h"
#include "veclib/veclib.h"

struct coef_table {                /* .... Integration coefficients .... */
  double gamma                  ;  /* implicit coefficient               */
  double alpha [_MAX_TORDER]    ;  /* multi-step coefficients for u(t)   */
  double beta  [_MAX_TORDER]    ;  /* multi-step coefficients for f(t)   */
};                                 /* .................................. */

#define Reshuffle(U,N)\
  {  int J = N-1; Element *t  = U[J];\
     for (q = J; q; q--) U[q] = U[q-1];\
     *U = t;  }

/* Adams-Bashforth schemes up to third order */

static struct coef_table AB_coef[] = {
  { 1., { 1.,   0.,  0. }, {   1.   ,    0. ,     0. }},
  { 1., { 1.,   0.,  0. }, {  3./2. , -1./2.,     0. }},
  { 1., { 1.,   0.,  0. }, { 23./12., -4./3., 5./12. }}
};

/* Stiffly-stable schemes up to third order */

static struct coef_table SS_coef[] = {
  {  1.,    { 1.,  0.,    0.   },  { 1.,  0.,  0. }},
  {  3./2., { 2., -1./2., 0.   },  { 2., -1.,  0. }},
  { 11./6., { 3., -3./2., 1./3.},  { 3., -3.,  1. }},
};

static struct coef_table *Icoef  = SS_coef;

/* ------------------------------------------------------------------------ *
 * Integrate() - ODE Integrator for U' = F(U,t)                             *
 *                                                                          *
 * This function implements a generic integration scheme for multi-level    *
 * time stepping.  The integration is of the form...                        *
 *                                                                          *
 *                   Je-1                           Je-1                    *
 * gamma U[n+1] = SUM (alpha[q] * U[n-q]) + dt * SUM (beta[q] * F[n-q])     *
 *                   q=0                            q=0                     *
 *                                                                          *
 * The coefficients alpha, beta, and gamma make the scheme accurate up to   *
 * order "J" in time, and can be either based on an Adams-Bashforth scheme  *
 * or a Stiffly-Stable scheme.                                              *
 *                                                                          *
 * The routine performs reshuffling of "Us" and "Fs" assuming they are of   *
 * the same dimension as the global parameter TORDER.                       *
 * ------------------------------------------------------------------------ */

double Alpha [_MAX_TORDER] = { 1., 0., 0. },
       Beta  [_MAX_TORDER] = { 1., 0., 0. },   
       Gamma               =   1.;

void Integrate (Field *U, Field *Us[], Field *Fs[], int Je, double dt)
{
  int      ntotz = U->nr * U->ns * U->nz * Field_count(U),
           Jmax  = iparam("TORDER");
  double   BetaDt [_MAX_TORDER];
  register int q;

  /* Clamp the value of Je to the range (1, Jmax), and get the coefficients *
   * for the integration from the table of stored values.  Also, set gamma  *
   * to the proper value for time stepping with order Je.                   */
  
  Je     = CLAMP(Je, 1, Jmax);
  Gamma  = Icoef[Je-1].gamma;
  for (q = 0; q < Je; q++) {
    Alpha  [q] = Icoef[Je-1].alpha[q];
    BetaDt [q] = (Beta[q] = Icoef[Je-1].beta[q]) * dt;
  }


  /* Prepare the solution array */

  dcopy (ntotz,         *U->base, 1, *Us[0]->base, 1);
  dscal (ntotz, *Alpha, *U->base, 1);


  /* Do the integration */
  
  for (q = 1; q < Je; q++)
    daxpy (ntotz, Alpha [q], *Us[q]->base, 1, *U->base, 1);
  for (q = 0; q < Je; q++)
    daxpy (ntotz, BetaDt[q], *Fs[q]->base, 1, *U->base, 1);
  

  /* Reshuffle the multi-step arrays */
  
  Reshuffle (Us, Jmax);
  Reshuffle (Fs, Jmax);

  return;
}

/* ------------------------------------------------------------------------ *
 * Because the integration coefficients are maintained statically within    *
 * this file, the following functions allow access to them from outside.    *
 * ------------------------------------------------------------------------ */

void set_order (int Je)
{
  const int Jmax = iparam("TORDER");

  Je    = CLAMP (Je, 1, Jmax) - 1;
  Gamma = Icoef[Je].gamma;

  memcpy (Alpha, Icoef[Je].alpha, _MAX_TORDER * sizeof(double));
  memcpy (Beta , Icoef[Je].beta , _MAX_TORDER * sizeof(double));

  return;
}

void get_alpha (double *a)
{
  memcpy (a, Alpha, _MAX_TORDER * sizeof(double));
  return;
}

void get_beta (double *b)
{
  memcpy (b, Beta, _MAX_TORDER * sizeof(double));
  return;
}

double get_gamma (int Je) 
{ 
  if (Je == 0)
    return Gamma;
  else
    return Icoef[Je-1].gamma; 
}

void set_Itype (int type)
{
  switch (type) {
  case StifflyStable:
    Icoef = SS_coef;
    break;
  case AdamsBashforth:
    Icoef = AB_coef;
    break;
  default:
    Prism_error("Prism: invalid integration type in set_Itype\n");
    break;
  }

  return;
}

/* ------------------------------------------------------------------------- *
 * estimate_CFL() - Estimate the CFL number                                  *
 *                                                                           *
 * This function computes an estimate of the CFL number based on the mesh    *
 * spacing and local velocity. This is a fairly expensive operation so it    *
 * should only be done every once in a while.                                *
 *                                                                           *
 * The CFL number is computed as follows:                                    *
 *                                                                           *
 *     CFL_x   = min { |u| dt/dx }                                           *
 *     CFL_y   = min { |v| dt/dy }                                           *
 *     CFL_z   = min { |w| dt/dz }                                           *
 *                                                                           *
 * where                                                                     *
 *                                                                           *
 *        dx   =  (dx/dr * dr + dx/ds * ds)                                  *
 *        dy   =  (dy/dr * dr + dy/ds * ds)                                  *
 *                                                                           *
 * ------------------------------------------------------------------------- */

#define SIGMA    0.9       /* safety factor                  */
#define CFL_min  0.25      /* smallest practical CFL number  */
#define CFL_max  0.70      /* stability limit for the scheme */

static  double   dzmin();  /* returns min mesh spacing */

void estimate_CFL (Domain *omega)
{
  Field *U = omega->U;
  Field *V = omega->V;
  Field *W = omega->W;  

  const int nrns = FIELD_NR(U)*FIELD_NS(U);
  const int nz   = U->nz;
  const int nel  = Field_count(U);

  double  dz, dr, ds, dt, CFL_x, CFL_y, CFL_z, CFL_dt, dt_max;
  int i, k;

  tempVector (tmp, nrns);

  dr    = dzmin(U->nr);
  ds    = dzmin(U->ns);
  dt    = dparam("DT");
  dz    = scalar("2*PI/(BETA*NZ)");
  CFL_x = CFL_y = CFL_z = 0.;

  for (k = 0; k < nel; k++) {
    for (i = 0; i < nrns; i++)   tmp[i]  = dr * fabs(U[k].xr[0][i]);
    if (U[k].xs)
      for (i = 0; i < nrns; i++) tmp[i] += ds * fabs(U[k].xs[0][i]);

    dvdiv (nrns, *U[k].base, 1, tmp, 1, tmp, 1);
    i     = idamax (nrns, tmp, 1);
    CFL_x = MAX (CFL_x, fabs(tmp[i]));

    for (i = 0; i < nrns; i++)   tmp[i]  = ds * fabs(U[k].ys[0][i]);
    if  (U[k].yr)
      for (i = 0; i < nrns; i++) tmp[i] += dr * fabs(U[k].yr[0][i]);

    dvdiv(nrns, *V[k].base, 1, tmp, 1, tmp, 1);
    i     = idamax (nrns, tmp, 1);
    CFL_y = MAX (CFL_y, fabs(tmp[i]));
  }
  
  CFL_dt = MAX (CFL_x, CFL_y);
  CFL_dt = MAX (CFL_z, CFL_dt);
  dt_max = SIGMA * CFL_max / CFL_dt;
  
  Prism_log(0,"\tCFL = %g, (dt) max = %g, dt = %g [%.0f%%]\n",
	    CFL_dt * dt, dt_max, dt, 100.*dt/dt_max);

  freeVector (tmp);
}

static double dzmin(int n)
{
  double *z;
  getops (n, &z, 0, 0, 0);
  return z[1]-z[0];
}
