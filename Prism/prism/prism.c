/*
 *                            P  R  I  S  M 
 *
 *                   Spectral Element / Fourier Solver
 *
 *      A Simulator for the Incompressible Navier-Stokes Equations
 *
 * The following program solves the 3D incompressible Navier-Stokes equations
 * using a combination of two-dimensional spectral elements and Fourier 
 * modes in the homogeneous z-direction.  The Navier-Stokes equations are 
 * treated as a set of evolution equations for the velocity field U:
 *
 *         d U        1          1
 * (1)     ---   =  - - grad p + - grad (grad U) - N(U) + F
 *         d t       rho         Re
 *
 * where U is the velocity field, rho (= 1) is the density, p is the static 
 * pressure, F is an applied force, and Re is the Reynolds number based on 
 * appropriate length and velocity scales.  
 * 
 * The simulation is constrolled by setting parameters in the input file. 
 * This includes both physical parameters related to the flow and numerical 
 * parameters that define the resolution and accuracy of the solution.  The
 * most important are:
 *
 *   Type        Name     Def   Description
 * --------     ------   -----  ------------------------------------------
 * [float]	Re		Reynolds number
 * [float]	TIME		Integration time interval
 *
 * [float]	DT       0.001	Time step
 * [integer]	NORDER	 5	Order + 1 of the (x,y)-basis functions
 * [integer]	NZ	 1	Number of frames (2 x modes) for 3D flows
 * [integer]    EQTYPE	 0	Form of the advection terms (see below)
 * [float]	TOLREL	 1.e-8	Error tolerance for all systems A u = F
 * [integer] 	MWORDS	 0	Maximum memory for matrices, 0 = no limit
 * 
 * Eq. (1) is integrated using a three-step splitting scheme.  The pressure
 * and diffusion equations are solved implicitly, while the nonlinear term is
 * integrated explicitly.  The program operates in rotational mode, skew-
 * symmetric mode, or Stokes mode:
 *
 *    Rotational 	EQTYPE = 0    N(U) =  U x grad x U	
 *    Skew-Symmetric 	EQTYPE = 1    N(U) = (U.grad U + grad.(UU)) * (1/2)
 *    Stokes     	EQTYPE = 2    N(U) =  0
 *
 * All calculations are double precision, and by default all systems are
 * solved directly.  If memory is limited, use the command line switch -i 
 * to use the iterative (PCG) solver for the viscous diffusion equations.
 *
 * Program development by:
 *
 * Ron Henderson
 * Department of Mechanical and Aerospace Engineering
 * Princeton University
 * Princeton, NJ 08540
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prism/prism.h"
#include "prism/config.h"
#include "veclib/veclib.h"

static char *prog  = "prism";
static char *usage = "[options] session[.rea]";

/* ------------------------------------------------------------------------- */

int main (int argc, char *argv[])
{ 
  double t0;
  double tN;
  Domain *domain;
  char   *session;

  /* Initialize */

  Prism_init(&argc,&argv);
  parse_args( argc, argv);
  session = argv[argc-1];

  /* Create the computational domain and get integration limits */

  domain = Domain_alloc(session);
  t0     = dparam("TIME_0");
  tN     = dparam("TIME_N");

  /* This is the simplest case possible: we just want to integrate the  *
   * Navier-Stokes equations from time t0 -> tN.  The "state" of the    *
   * fluid is defined by the initial conditions, so now we want to see  *
   * what happens as it moves to the new state at t = tN.               */

  Navier_Stokes (domain, t0, tN);
  PostProcess   (domain);

  /* The End */

  Domain_free(domain);
  Prism_exit();
}


