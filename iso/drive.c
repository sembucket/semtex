/*****************************************************************************
 * Iso:
 * ----
 * DNS of isotropic turbulence by a Fourier--Galerkin pseudo-spectral method.
 *
 * Copyright (C) 1992-1999 Hugh Blackburn
 *
 * Usage:
 * ------
 * iso [options] session
 * options:
 * -h       ... print this message
 * -c       ... turn off checkpointing         [Default: on]
 * -d <num> ... set time step                  [Default: 0.01]
 * -f <num> ... set interval for field dumps   [Default: 100]
 * -H <num> ... set interval for history dumps [Default: 10]
 * -k <num> ... set kinematic viscosity        [Default: 1]
 * -n <num> ... set number of timesteps        [Default: 100]
 * -o <num> ... set timestepping order         [Default: 2]
 *
 * Method:
 * -------
 * Time evolution of Fourier coefficients of flow in box with sides of length
 * 2PI and periodic boundary conditions.  Number of modes in each direction
 * must be the same integer power of two.
 *
 * Computation of the non-linear terms is fully de-aliased by the method
 * of isotropically-truncated convolution sums in conjunction with phase
 * shifts, as described by Orszag [1].
 *
 * Time advancement is explicit, Adams-Bashforth 2 on the nonlinear
 * terms, with an analytic integrating-factor treatment of viscous terms,
 * developed by Rogallo [2].  Initial timestep is Euler.
 *
 * A general description of the Fourier treatment and elimination of
 * pressure by projection in Fourier-space can be found in ch. 4 of
 * Lesieur's book [3].
 *
 * Storage requirements:
 * ---------------------
 * In terms of main storage, the programme carries 14*N*N*N words, with
 * additional incidental requirements for storage of one-dimensional vectors,
 * pointers and scalar variables.
 *
 * Main storage consists of the Fourier coefficients of the velocities
 * (3*N*N*N), the nonlinear terms at the current timestep (G) and the last
 * timestep (G_old) (6*N*N*N), and workspace (5*N*N*N).  The most obvious
 * ways to save storage would be to carry only two velocity components, and
 * to keep G_old on disk: this would produce a saving of 4*N*N*N.  With much
 * hard work, an additional N*N*N of workspace could probably be saved.
 * For large simulations, out-of-core storage will be needed.
 *
 * Workload:
 * ---------
 * The bulk of CPU time is spent in 3D real--complex FFTs.  For the fully-
 * dealiased computations, there are 21 FFTs per timestep.  For a solution
 * with aliasing, this could be reduced to 12 (storage requirements would
 * also be reduced to 10*N*N*N words --- 9 for velocity fields and nonlinear
 * terms, 1 for workspace).
 *
 * Files:
 * ------
 * A binary field file that gives initial conditions is required.
 * All field-type files are written as Fourier coefficients.
 *
 * Naming conventions for files (suffixes are appended to session file root):
 *   .rst   Restart/IC file
 *   .fld   Field      file     
 *   .chk   Checkpoint file 
 *
 * References:
 * -----------
 * [1]  S.A. Orszag, 1971, 'Numerical simulation of incompressible flows
 *        within simple boundaries. 1. Galerkin (spectral) representations',
 *        Stud. Appl. Math., VLN4, Dec., 293--327.
 * [2]  R.S. Rogallo, 1981, 'Numerical experiments in homogeneous
 *        turbulence', NASA TM 81315, Sept.
 * [3]  M. Lesieur, 1990,  Turbulence in Fluids, 2nd ed, Kluwer Academic.
 * [4]  C. Canuto, M.Y. Hussaini, A. Quarteroni & T.A. Zang, 1988,  Spectral
 *        Methods in Fluid Dynamics, Springer.
 * [5]  M.E. Brachet, D. I. Meiron, S.A. Orszag, B.G. Nickel, R.H. Morf
 *        & U. Frisch, 1983.  'Small-scale structure of the Taylor--Green
 *        vortex', JFM, V130, 411--452.
 *
 * Program development by:
 * -----------------------
 * Hugh Blackburn
 * Department of Mechanical Engineering
 * Monash University
 * Clayton VIC 3168
 * Australia
 * hmb@artemis.eng.monash.edu.au
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

static const char prog[] = "iso";
static void Default (Param*);
static void getargs (int, char**, Param*);

int N, K, FourKon3;		/* -- Global grid size variables. */

#define ROLL(u, n) { if ((n) > 1) { CVF _tmp = (u)[(n)-1]; int _i; \
  for (_i = (n)-1; _i; _i--) (u)[_i]=(u)[_i-1]; (u)[0] = _tmp; }}


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver routine for simulation code.
 * ------------------------------------------------------------------------- */
{
  CVF     U, work1;
  CVF*    G = malloc (sizeof (CVF*));
  CF      work2, work3;
  complex *Wtab, *Stab;
  Param*  I = calloc (1, sizeof (Param));
  real    setKinvis;

  /* -- Set up. */

  Default (I);
  getargs (argc, argv, I);
  setKinvis = I -> kinvis;
  startup (I);

  printf ("%s: Fourier-spectral Navier-Stokes simulation\n", prog);
  printf ("Copyright (C) 1992-1999 Hugh Blackburn\n\n");

  N = I -> ngrid; K = N / 2; FourKon3 = (4 * K) / 3;

  allocate   (&U, G, I -> norder, &work1, &work2, &work3, &Wtab, &Stab);

  preFFT     (Wtab, K);
  preShift   (Stab, N);

  restart    (U, I);
  I -> kinvis = setKinvis;

  truncateVF (U);
  projectVF  (U, work2);

  printParam (stdout, I);
  analyze    (U, I, Wtab);

  /* -- Time-stepping loop. */

  while (I -> step < I -> nstep) {

    zeroVF    (G[0]);
    nonlinear (U, G[0], work1, work2, work3, Wtab, Stab);
    projectVF (G[0], work2);
    integrate (U, G, I);
    ROLL      (G, I -> norder);

    I -> step ++;
    I -> time += I -> dt;
    
    if (!(I -> step % 100)) projectVF (U, work2);

    analyze (U, I, Wtab);
    dump    (U, I);
  }

  cleanup (I);

  return EXIT_SUCCESS;
}


static void Default (Param* I)
/* ------------------------------------------------------------------------- *
 * Set run-time defaults.
 * ------------------------------------------------------------------------- */
{
  I -> fld_dmp = 0;
  I -> his_dmp = 0;
  I -> io_fld  = 100;
  I -> io_his  = 10;
  I -> chkpnt  = TRUE;
  I -> ngrid   = 16;
  I -> nstep   = I -> io_fld;
  I -> norder  = 2;
  I -> step    = 0;
  I -> dt      = 0.01;
  I -> time    = 0.0;
  I -> kinvis  = 1.0;
}


static void getargs (int    argc,
		     char** argv,
		     Param* I   )
/* ------------------------------------------------------------------------- *
 * Process command-line arguments.
 * ------------------------------------------------------------------------- */
{
  const char usage[] =
    "usage: %s [options] session\n"
    "options:\n"
    "-h       ... print this message\n"
    "-c       ... turn off checkpointing         [Default: on]\n"
    "-d <num> ... set time step                  [Default: 0.01]\n"
    "-f <num> ... set interval for field dumps   [Default: 100]\n"
    "-H <num> ... set interval for history dumps [Default: 10]\n"
    "-k <num> ... set kinematic viscosity        [Default: 1]\n"
    "-n <num> ... set number of timesteps        [Default: 100]\n"
    "-o <num> ... set timestepping order         [Default: 2]\n";
  char c, err[STR_MAX];
 
  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fprintf (stdout, usage, prog);
      exit (EXIT_SUCCESS);
      break;
    case 'c':
      I -> chkpnt = FALSE;
      break;
    case 'd':
      if   (*++argv[0]) I -> dt     = atof (*argv);
      else { --argc;    I -> dt     = atof (*++argv); }
      break;
    case 'f':
      if   (*++argv[0]) I -> io_fld = atoi (*argv);
      else { --argc;    I -> io_fld = atoi (*++argv); }
      break;
    case 'H':
      if   (*++argv[0]) I -> io_his = atoi (*argv);
      else { --argc;    I -> io_his = atoi (*++argv); }
      break;
    case 'k':
      if   (*++argv[0]) I -> kinvis = atof (*argv);
      else { --argc;    I -> kinvis = atof (*++argv); }
      break;
    case 'n':
      if   (*++argv[0]) I -> nstep  = atoi (*argv);
      else { --argc;    I -> nstep  = atoi (*++argv); }
      break;
    case 'o':
      if   (*++argv[0]) I -> norder = atoi (*argv);
      else { --argc;    I -> norder = atoi (*++argv); }
      break;
    default:
      break;
    }

  if (I->norder > 3) message (prog, "maximum timestepping order is 3", ERROR);

  if   (argc != 1)  message (prog, "no session definition file", ERROR);
  else              I -> session = *argv;
}
