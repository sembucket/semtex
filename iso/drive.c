/*****************************************************************************
 * Iso:
 * ----
 * DNS of isotropic turbulence by a Fourier--Galerkin pseudo-spectral method.
 *
 * Usage:
 * ------
 * iso [options] session
 * options:
 *  -h   ... print this message
 *  -chk ... checkpoint field dumps
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
 * Two files are needed; an ASCII session file, which gives run
 * parameters, and a binary field file which gives initial conditions.
 * All field-type files are written as Fourier coefficients.
 *
 * Naming conventions for files (suffixes are appended to session file root):
 *   .rst   Restart/IC file
 *   .fld   Field      file     
 *   .chk   Checkpoint file 
 *
 * Example session file:            (order is significant)
 *   Test run                TITLE  (tag string, otherwise unused)
 *   64                      NGRID  (spatial resolution)
 *   0.0001                  DT     (time step)
 *   1000                    IO_FLD (steps between fields or checkpoints)
 *   20000                   STEPS  (total number of steps to run)
 *   500.0                   RE     (Reynolds number: 1/(kinematic viscosity))
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


#define  SWAP(a, b) {CVF G_temp = (a); (a) = (b); (b) = G_temp;}

static void getargs (int, char**, char**, int*);
static char prog[] = "iso";


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver routine for simulation code.
 * ------------------------------------------------------------------------- */
{
  CVF      U, G, G_old, work;
  CF       F, F_;
  int*     Dim  = ivector (1, 3);
  complex* Wtab;
  complex* Stab;
  Param*   runInfo = calloc (1, sizeof (Param));
  char*    session;
  FILE*    fp;
  int      chkpoint = FALSE;

  /* -- Set up. */

  getargs (argc, argv, &session, &chkpoint);

  fp = efopen (session, "r");
  startup     (fp, runInfo, session, chkpoint);
  fclose      (fp);
  printParam  (stdout, runInfo, "$RCSfile$", "$Revision$");

  Dim[1] = (Dim[2] = runInfo -> modes);
  Dim[3] =  Dim[1] / 2;

  allocate   (&U, &G, &G_old, &work, &F, &F_, &Wtab, &Stab, Dim);
  initialize (U, runInfo, Dim);

  preFFT   (Wtab, Dim[3]);
  preShift (Stab, Dim[1]);

  analyze  (U, runInfo, Wtab, Dim);

  /* -- Main time-stepping loop. */

  while (runInfo -> step < runInfo -> stepMax) {

    nonlinear (U, G, F, F_, work, Wtab, Stab, Dim);
    integrate (U, G, G_old, runInfo, Dim);
    SWAP      (G, G_old);

    runInfo -> step ++;
    runInfo -> time += runInfo -> dt;

    analyze (U, runInfo, Wtab,     Dim);
    dump    (U, runInfo, chkpoint, Dim);
  }

  cleanup (runInfo, chkpoint);

  return EXIT_SUCCESS;
}


static void getargs (int argc, char **argv, char **session, int *chkpnt)
/* ------------------------------------------------------------------------- *
 * Process command-line arguments.
 * ------------------------------------------------------------------------- */
{
  char c, routine[] = "getargs";
  char usage[]   =
    "Usage: %s [options] session\n"
    "  [options]:\n"
    "  -h   ... print this message\n"
    "  -chk ... checkpoint field dumps\n";
 
  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fprintf (stdout, usage, prog);
      exit (EXIT_SUCCESS);
      break;
    case 'c':
      if (strstr ("chk", *argv)) {
	*chkpnt = TRUE;
      } else {
	fprintf (stdout, usage, prog);
	exit (EXIT_FAILURE);	  
      }
      break;
    default:
      break;
    }

  if   (argc != 1)  message (routine, "no session definition file", ERROR);
  else             *session = *argv;
}
