/*
 * Driver for the Navier-Stokes Solver
 *
 * This is the top-level integration function for the Navier-Stokes simulation.
 * It integrates the fluid system defined in "Omega" from its current state to
 * a new state at time t = t0 + dt*NSTEPS.  
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "veclib/veclib.h"
#include "prism/prism.h"

/* Hooks */

void (*user_advect)(Domain *);



static int    direct, flowrt;     /* Flags */

static int    Je, pid;            /* Externals */
static double dt, Re, *P_save;
static double tolu, tolp;
static char   rcsid[] = "$Revision$";

/* ------------------------------------------------------------------------- */

int Navier_Stokes (Domain *omega, double t0, double tN)
{ 
  int          step, nsteps;      /* Discrete time (step)       */
  double       time;              /* Continuous time            */
  ACTION       WaveProp;          /* Convective form            */
  Field        *V  [DIM],         /* Velocity field             */
              **Vs [DIM],         /* Multi-step storage         */
              **Vf [DIM];         /* Forcing functions          */
  Bedge        *Vbc[DIM], *Pbc;   /* Boundary conditions        */
  BSystem      *Mv,  *Mp;         /* Global matrix systems      */
  register int  i;

  Scalar       P;                 /* Pressure                   */


  dt       = dparam ("DT");                    /* Time Stepping */
  nsteps   = (tN - t0 + dt/2.) / dt;
  dparam_set ("TIME", time = t0);
  iparam_set ("STEP", step =  0);

  tolu     = dparam ("TOLU");          /* Simulation parameters */
  tolp     = dparam ("TOLP");
  Re       = dparam ("Re");
  Je       = iparam ("TORDER");
  WaveProp = (ACTION) iparam ("EQTYPE");
  pid      = option ("procid");
  direct   = option ("direct");
  flowrt   = option ("flowrate");

  Mv       = omega->Velocity;                 /* Matrix Systems */
  Mp       = omega->Pressure;

  P        = omega->P;            /* Pressure & Velocity fields */
  Pbc      = omega->Pbc;
  V  [0]   = omega->U  ; V  [1] = omega->V  ;
  Vbc[0]   = omega->Ubc; Vbc[1] = omega->Vbc;
  Vf [0]   = omega->Uf ; Vf [1] = omega->Vf ;
  Vs [0]   = omega->Us ; Vs [1] = omega->Vs ;
# if DIM == 3
  V  [2]   = omega->W  ;
  Vbc[2]   = omega->Wbc;
  Vf [2]   = omega->Wf ;
  Vs [2]   = omega->Ws ;

  if (!direct && !P_save)           /* Extra storage for pressure */
    P_save = (double*) 
      calloc (V[0]->nr * V[0]->ns * V[0]->nz * Mv->elements, sizeof(double));
#endif

  /* If the flowrate option is turned on, compute the Green's functions  *
   * for the Flowrate() routine and store them in Omega.                 */

  ROOTONLY
    if (option("flowrate") && !omega->G) flowrate_init (omega);

  Prism_log(0,"\nMarch from t=%g to %g using dt=%g [%d steps]\n", 
	    t0, tN, dt, nsteps);

  Vorticity (omega);                 /* Initialize vorticity */
  Analyzer  (omega, time, step);     /* Analyzer startup     */

  /* ======================================================= */
  /*                                                         */
  /*              MAIN   INTEGRATION   LOOP                  */
  /*                                                         */
  /* ======================================================= */

  while (step < nsteps) {

    MakeF(omega, Prep, step += 1, time);

    /* ..........    Advection step    .......... */

    MakeF(omega, WaveProp, step, time);
    for  (i = 0; i < DIM; i++) 
      Integrate (V[i], Vs[i], Vf[i], step, dt);

    /* ..........     Pressure step     .......... */

    MakeF      (omega, Pressure, step, time);
    MultiSolve (omega, &P, Vf, &Pbc, Mp, 1, step);

    /* ..........  Viscous correction   .......... */

    MakeF      (omega, Viscous, step, time);
    MultiSolve (omega, V, Vf, Vbc, Mv, DIM, step);
    
    MakeF      (omega, Post, step, time += dt);
    Analyzer   (omega, time, step);
  }

  return 0;   /* Normal exit */
}

/* ------------------------------------------------------------------------- *
 * MakeF() - Create a Forcing Function                                       *
 *                                                                           *
 * This function creates a generic right-hand-side vector field.  The input  *
 * is a pointer to the current domain and a requested "action" which deter-  *
 * mines which fields should be changed.  The possible choices are:          *
 *                                                                           *
 * Prep         -                                                            *
 *                                                                           *
 * Rotational   - Compute the nonlinear terms in U x (grad x U) form. For    *
 *                this formulation, P becomes total pressure.                *
 *                                                                           *
 * SkewSymm.    - Compute the nonlinear terms in skew-symmetric form.        *
 *                                                                           *
 * Stokes       - Add only the forcing function, i.e., N(U) = 0.             *
 *                                                                           *
 * Pressure     - Compute div U', where U' is the first intermediate veloc-  *
 *                ity field (the non-linear part).  This is the RHS of the   *
 *                pressure Poisson equation to correct the divergence.       *
 *                                                                           *
 * Viscous      - Compute U'', the second intermediate velocity field.  This *
 *                the first velocity field corrected for divergence errors.  *
 *                                                                           *
 * Post         - Do POST calculations to wrap up the current time step. If  *
 *                "flowrate" is turned on, this corrects the velocity field  *
 *                for the proper flowrate. Then it computes the vorticity    *
 *                of the current velocity field for analysis and for use in  *
 *                the next time step.                                        *
 * ------------------------------------------------------------------------- */

void MakeF (Domain *omega, ACTION act, int step, double time)
{
  Field   *Uf   = *omega->Uf;
  Field   *Vf   = *omega->Vf;
  Field   *Wf   = *omega->Wf;
  Field   *Pf   = *omega->Uf;
  BSystem *M    =  omega->Pressure;

  const int nz   =  Uf->nz;
  const int ntot =  Uf->nr * Uf->ns * M->elements;
          
  switch (act) {

  case Prep: 
    iparam_set ("STEP", step);      /* New time step            */
    set_order  (CLAMP(step,1,Je));  /* Update integration table */
    break;

  case Rotational:
    VxOmega (omega);  goto AddForcing; 
    
  case SkewSymmetric:
    SkewSymm(omega);  goto AddForcing; 

  case Stokes:
    StokesBC (omega); goto AddForcing;

  case Convective:
    if (user_advect) {
      (*user_advect)(omega); goto AddForcing;
    } else
      Prism_error("user_advect: undefined!");

  AddForcing: {
    double f;

    /* Standard Forcing (constant) */

    ROOTONLY {
      if ((f=dparam("FFX")) != 0.)
	dsadd(ntot, f, *Uf->base, 1, *Uf->base, 1);
      if ((f=dparam("FFY")) != 0.)
	dsadd(ntot, f, *Vf->base, 1, *Vf->base, 1);
      if ((f=dparam("FFZ")) != 0.)
	dsadd(ntot, f, *Wf->base, 1, *Wf->base, 1);
    }

    SetPBCs(omega);
    break;
  }
    
  case Pressure: {
    Field    *U    = omega->U,
             *V    = omega->V;
    double   dtinv = 1./dt;
    register int i, m;
    
  /* Apply a DSSUM to U to finish the advection step */

# if DIM == 3
#
    Field    *W    = omega->W;

    Field_davg_3D (U, omega->Velocity);
    Field_davg_3D (V, omega->Velocity);
    Field_davg_3D (W, omega->Velocity);

    Field_grad_3D (U, Uf, 0, *omega->Us[0]->base, *omega->Vs[0]->base);
    Field_grad_3D (V, 0, Vf, *omega->Us[0]->base, *omega->Vs[0]->base);
    Field_gradz   (W, Wf);

    for (m = 0; m < nz; m = m || pid ? m+1 : m+2)
      for (i = 0; i < ntot; i++)
	(*Pf->base)[i + m*ntot] = dtinv * 
	  ((*Uf->base)[i+m*ntot]+(*Vf->base)[i+m*ntot]+(*Wf->base)[i+m*ntot]);
    
#
# else
#
    Field_davg (U, omega->Velocity);
    Field_davg (V, omega->Velocity);

    Field_grad (U, Uf, 0);
    Field_grad (V, 0, Vf);

    dsvvpt (ntot, dtinv, *Uf->base, 1, *Vf->base, 1, *Pf->base, 1);
#
#endif
    if (!direct) dcopy (ntot*nz, P_save, 1, *omega->P->base, 1);
    break;
  }

  case Viscous: {
    Field    *U    = omega->U,
             *V    = omega->V,
             *W    = omega->W,
             *P    = omega->P;
    register int i, m;

    /*
    //  NB: Solve the "advection equation" 
    //
    //        u^^ - u^ = -dt * grad P,
    //
    //      including the DSSUM to get the correct u^^.  
    //
    //      DSSUM:  Don't leave home without it!
    */

#
# if DIM == 3
#
    if (!direct) dcopy (ntot*nz, *P->base, 1, P_save, 1);

    Field_grad_3D (P,  Uf, Vf, *omega->Us[0]->base, *omega->Vs[0]->base);
    Field_davg_3D (Uf, omega->Velocity);
    Field_davg_3D (Vf, omega->Velocity);
    Field_gradz   (P,  Wf);

    for (i = m = 0; m < nz; m = m || pid ? m+1 : m+2, i = m*ntot) {
      daxpy (ntot, -dt,   *Uf->base + i, 1,  *U->base + i, 1);
      daxpy (ntot, -dt,   *Vf->base + i, 1,  *V->base + i, 1);
      daxpy (ntot, -dt,   *Wf->base + i, 1,  *W->base + i, 1);

      dsmul (ntot, -Re/dt, *U->base + i, 1, *Uf->base + i, 1);
      dsmul (ntot, -Re/dt, *V->base + i, 1, *Vf->base + i, 1);
      dsmul (ntot, -Re/dt, *W->base + i, 1, *Wf->base + i, 1);
    }
#
#else
#
    Field_grad (P,  Uf, Vf);
    Field_davg (Uf, omega->Velocity);
    Field_davg (Vf, omega->Velocity);

    daxpy (ntot, -dt, *Uf->base, 1, *U->base, 1);
    daxpy (ntot, -dt, *Vf->base, 1, *V->base, 1);
    
    dsmul (ntot, -Re/dt, *U->base, 1, *Uf->base, 1);
    dsmul (ntot, -Re/dt, *V->base, 1, *Vf->base, 1);
#
#endif

    break;
  }
    
  case Post: 
    ROOTONLY if (flowrt) flowrate(omega);       /* Flowrate */

    if (!option("need.vorticity")) 
      Vorticity (omega);                        /* Vorticity */

    dparam_set ("TIME", time);                  /* Update time */
    break;

  default:
    Prism_error("MakeF: unknown type of ACTION\n");
    break;
  }
  
  return;
}

/* ------------------------------------------------------------------------- *
 * MultiSolve() - Wrapper for Solve                                          *
 *                                                                           *
 * This is a wrapper for the Helmholtz solver to compute the solution for a  *
 * number of fields with the same operator.  It also decides whether to use  *
 * the direct or iterative solver based on the type of field (velocity or    *
 * pressure) and the value of the parameter "direct".  There are three pos-  *
 * sible actions:                                                            *
 *                                                                           *
 *       direct = 0  mean pressure direct, everything else iterative         *
 *                1  pressure direct, velocity iterative                     *
 *                2  pressure direct, velocity direct if Je <= step          *
 * ------------------------------------------------------------------------- */

void MultiSolve
  (Domain *omega, Field *U[], Field **F[], Bedge *Ubc[], BSystem *B, 
                  int nfields,  int step)
{
  char      type  = (*U)->type;
  const int nz    = (*U)->nz;
  double    beta, lambda, wavek;
  register int k, n, mz;

  step   = MIN (step, Je);
  beta   = dparam("BETA");
  wavek  = (k = 0);

  if (type == 'u') {
    dparam_set("TOLREL", tolu);
    B->constant = lambda = get_gamma(step) * Re / dt;
  } else {
    dparam_set("TOLREL", tolp);
    B->constant = lambda = 0.;
  }

  /* ----------------------------------------------------- *
   *       2 - D   Velocity and Pressure Solver            *
   * ----------------------------------------------------- */
#
# if DIM == 2     
#
  B->constant = lambda;
  if (type == 'u' && (Je > step || direct < 2)) {
    for (n = 0; n < nfields; n++) 
      Solve_CG(U[n], *F[n], Ubc[n], B);
  } else
    for (n = 0; n < nfields; n++) 
      Solve   (U[n], *F[n], Ubc[n], B);
#
# else 
#
  /* ----------------------------------------------------- *
   *       3 - D   Velocity and Pressure Solver            *
   * ----------------------------------------------------- */

  if (direct < 2 || Je > step) {
# ifndef PCG
    if (type == 'p' && pid == 0) {                  /* one DIRECT solve */
      Frame_load (*U, B, k);
      Frame_set  (k, 2,  U[0],  *F[0]);
      Solve      (U[0], *F[0], Ubc[0], B);
      k += 2;
    }
# endif
    
    if (type != 'p' || direct < 1) {  /* ...then the rest are ITERATIVE */
      while (k < nz) {
	if (!(k & 1)) {
	  wavek     = .5 * beta * (k + pid*nz);
	  B->constant = lambda + wavek*wavek;
	}

	for (n = 0; n < nfields; n++) {
	  Frame_set(k, 2,  U[n],  *F[n]);
	  Solve_CG (U[n], *F[n], Ubc[n], B);
	}
	k = (k || pid) ? k+1 : k+2;
      }
    }  /* NOTE: at the exit of this loop, k = nz */
  }   
 

  /* --------  All DIRECT solves  -------- */

  while (k < nz) {

    if (!(k & 1)) {
      wavek     = .5 * beta * (k + pid*nz);
      B->constant = lambda + wavek*wavek;
      Frame_load (*U, B, (k+pid*nz)/2);
    }

    for (n = 0; n < nfields; n++) {
      Frame_set(k, 2,  U[n],  *F[n]);
      Solve    (U[n], *F[n], Ubc[n], B);
    }

    k = (k || pid) ? k+1 : k+2;
  }

  /* Set the imaginary part of mode 0 to zero */
  
  ROOTONLY {
    const int ntot = (*U)->nr * (*U)->ns * B->elements;
    for (n = 0; n < nfields; n++)
      dzero (ntot, *U[n]->base + ntot, 1);
  }
#
# endif
#
  return;
}
