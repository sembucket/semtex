/* 
 * Implementation of Domain
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "prism/prism.h"
#include "prism/domain.h"
#include "prism/measure.h"
#include "prism/hooks.h"
#include "speclib/speclib.h"

/* Hooks */

BSystem* (*user_build)(Field*, Bedge*, const char*, double);

/* Check all parameters defined in the input file */

static void check(void)
{
  int    n, nz;
  double val, dt;

  /* NB:  Only three dimensional parameters are required:  a velocity, 
  //      a length, and a density.   By default, these are assumed to
  //      have values of unity in arbitrary units (i.e. kg/m^3).  If
  //      you used something else do non-dimensionalize the mesh and
  //      boundary conditions, you can define VELOCITY, LENGTH, and
  //      DENSITY in the input file to get the appropriate scaling.
  //
  //      The Reynolds number, Re, is the only nondimensional parameter
  //      that appears in the Navier-Stokes equations.  Given Re and the
  //      basic velocity and length scales, we can generate dimensional
  //      values of the viscosity.  When computing forces that involve
  //      [dimensional] viscosity coefficients, be careful that you keep
  //      things properly nondimensional.
  //
  //      Example:  The [dimensional] shear stress is
  //
  //                tau_ij = DYNVIS * VELOCITY / LENGTH * 
  //                            ( d u_i / d x_j + d u_j / d x_i )
  //
  //                The drag coefficient for a body would be defined by
  //                integrating (n_x . tau_ij) over the body surface,
  //                and then normalizing it as
  //
  //                CD = ( Force ) / ( 1/2 * DENSITY * VELOCITY^2 )
  //
  //                Unless you do something funny, the force coefficient
  //                is simply 2 times the force evalulated directly from
  //                the simulation variables.
  */

  dparam_set ("DYNVIS",  scalar ("LENGTH * DENSITY * VELOCITY / Re"));
  dparam_set ("KINVIS",  scalar ("DYNVIS / DENSITY"));
	      
  /* The following section is only for unsteady problems */

  iparam_set ("TORDER", CLAMP(iparam("TORDER"), 1, _MAX_TORDER));

  if (!option("steady")) {
    if ((dt=dparam("DT")) == 0.) {
      if ((dt=dparam("DELT")) != 0.) 
	dparam_set("DT", dt);
      else
	Prism_error("Prism: no value specified for DT\n");
    }

    /* Round the integration time to a DT-multiple */

    if ((val=dparam("TIME")) > 0.)
      iparam_set("NSTEPS", n = (val+.5*dt)/dt);
    else
      val = dt * (n = iparam("NSTEPS"));

    dparam_set ("TIME",   val);    /* Interval length        */
    dparam_set ("TIME_0", 0.);     /* Initial time [default] */
    dparam_set ("TIME_N", val);    /* Final time             */

    /* Set default values */
    
    if (iparam("IO_FLD") == 0) iparam_set ("IO_FLD", n);
    if (iparam("IO_HIS") == 0) iparam_set ("IO_HIS", n / 1024);
    if (iparam("IO_MEA") == 0) iparam_set ("IO_MEA", n / 1024);

    /* Clamp the IO_{Class} parameters */

    iparam_set ("IO_FLD", CLAMP(iparam("IO_FLD"), 0, n));
    iparam_set ("IO_HIS", CLAMP(iparam("IO_HIS"), 0, n));
    iparam_set ("IO_MEA", CLAMP(iparam("IO_MEA"), 0, n));
    iparam_set ("IO_CFL", CLAMP(iparam("IO_CFL"), 0, n));
    
    /* Set the IO_{Class}_DT parameters */
    
    dparam_set ("IO_FLD_DT", dt * iparam("IO_FLD"));
    dparam_set ("IO_HIS_DT", dt * iparam("IO_HIS"));
    dparam_set ("IO_MEA_DT", dt * iparam("IO_MEA"));
    dparam_set ("IO_CFL_DT", dt * iparam("IO_CFL"));

    if (option("checkpt")) {        /* Default time for checkpoints is */
      n   /= 10 ;                   /* every 10% of the simulation.    */
      val /= 10.;
      iparam_set("IO_FLD",    MAX(n,1));
      dparam_set("IO_FLD_DT", MAX(val,dt));
    }
  }

  iparam_set ("DIM", DIM);

  if (dparam("LZ") > 0.) 
    dparam_set ("BETA", scalar("2*PI/LZ"));
  else if (dparam("BETA") > 0.)
    dparam_set ("LZ",   scalar("2*PI/BETA"));
  
  /* Convergence tolerances for the PCG-solver */

  val = dparam("TOLREL");
  if (dparam("TOLP") == 0.) dparam_set("TOLP", val);
  if (dparam("TOLU") == 0.) dparam_set("TOLU", val);

  /* Limiting memory usage ? */

  if (n = iparam("MWORDS")) {
    if (! iparam("MWORDSU")) iparam_set("MWORDSU", n);
    if (! iparam("MWORDSP")) iparam_set("MWORDSP", n);
  }
}


/* Set up matrix systems for the run */

static BSystem *buildops(Field *U, Bedge *Ubc, const char *name, double lambda)
{
  double  beta;
  BSystem *B;
  int     k, mz, nz, direct, pid, nprocs;

  nz     = U->nz;
  beta   = dparam("BETA");
  direct = option("direct");
  pid    = option("procid");
  nprocs = option("nprocs");

  B = Matrix_alloc (U, Ubc, name);

  /* Do a couple checks to see which matrices are needed for this field.  *
   * The three solve levels (direct = [0,1,2]) require (0) no velocity    *
   * matrices, 1 pressure matrix; (1) no velocity matrices, nz/2 pressure *
   * matrices; and (2) nz/2 velocity and pressure matrices.               */

  if (lambda != 0. && direct < 2) {
    Prism_log(0,"none");
    return B;
  }
# ifdef PCG
  else if (direct < 1) {
    Prism_log(0,"none");
    return B;
  }
# else
  else if (direct < 1) nz = (pid ? 0 : 1);
#endif  

  if (nprocs > 1 || nz > 32) 
    Prism_log(0,"...please wait...");

  for (k = 0; k < nz; k += 2) {

    mz = (k + pid*nz) >> 1;                    /* wave number in z   */
    B->constant = lambda + pow(beta*mz,2.);    /* Helmholtz constant */

    /* Try to load the matrix from a file.  If it doesn't exist then *
     * compute a new set and save it.  Handled by the frame routines */

    if (Frame_init (U, B, mz, name)) {

      Matrix_build (U, B);
      Frame_save   (U, B, mz);

      /* For small simulations echo a '*' for a banded system, a '.' for *
       * non-banded systems, and an 'f' if the matrix was succesfully    *
       * loaded from a file.  No echo for large simulations.             */

      if (nprocs > 1 || nz > 32) 
	continue;
      else
	Prism_log(0,((B->bandwidth > 0) ? "*" : "."));
    } else
      Prism_log(0,"f");
  }
  
#ifdef PARALLEL
  comm_sync();
#endif

  /* For 2-D simulations, the matrices always reside in memory, i.e., no  *
   * matrices are swapped from disk.  Since there is no call to LoadFrame *
   * in the main integration loop, the frames need to be loaded here once *
   * and for all.                                                         */

  if (DIM == 2) Frame_load (U, B, 0);

  return B;
}


/* ------------------------------------------------------------------------- */


/* Create a new computational domain */

Domain *Domain_alloc (const char *session)
{
  int    Je, n;
  double Re, dt;
  char   fname[FILENAME_MAX];
  FILE   *fp;

  /* Create the new domain and open the input file */

  Domain *domain = (Domain*) calloc(1,sizeof(Domain));
  domain->name   = strdup(session);

  sprintf (fname, "%s.rea", session);
  if (!(fp=fopen(fname,"r")))
    Prism_error("Domain_alloc: can't open the input file: %s\n", session);
  else
    Prism_log(0,"Reading problem specification from %s\n\n", fname);

  /* Read the input parameters */

  ReadParams(fp);
  check();

  domain->U = ReadMesh(fp);           FIELD_TYPE(domain->U) = 'u';
  domain->V = Field_dup (domain->U);  FIELD_TYPE(domain->V) = 'v';
#if DIM == 3
  domain->W = Field_dup (domain->U);  FIELD_TYPE(domain->W) = 'w';
#endif
  domain->P = Field_dup (domain->U);  FIELD_TYPE(domain->P) = 'p';

  /* Read the boundary conditions */

  domain->Ubc = ReadBCs(fp, 0, domain->U);
  domain->Vbc = ReadBCs(fp, 1, domain->V);
# if DIM == 3
  domain->Wbc = ReadBCs(fp, 2, domain->W);
# endif
  domain->Pbc = domain->Ubc ? BuildPBCs(domain->P, domain->Ubc) : NULL;
  
  /* Build the multistep arrays */

  if (!option("prep")) {
    Je = iparam("TORDER");
    for (n = 0; n < Je; n++) {
      (domain->Us)[n] = Field_dup(domain->U);
      (domain->Uf)[n] = Field_dup(domain->U);

      (domain->Vs)[n] = Field_dup(domain->V);
      (domain->Vf)[n] = Field_dup(domain->V);
#if DIM == 3
      (domain->Ws)[n] = Field_dup(domain->W);
      (domain->Wf)[n] = Field_dup(domain->W);
#endif
    }
  }

  /* Initialize the solution, drive force, history points, etc. */

  Summary     ();
  ReadICs     (fp, DIM, domain->U, domain->V, domain->W);
  ReadDF      (fp, DIM);
  ReadHisData (fp, domain);

  ROOTONLY {
    if (option("verbose") > 1) {
      Prism_options();
      Prism_params ();
    }
  }

  fclose (fp);

  /* Initialize output frequency parameters */

  domain->step.cfl     = iparam("IO_CFL");
  domain->step.history = iparam("IO_HIS");
  domain->step.measure = iparam("IO_MEA");
  domain->step.field   = iparam("IO_FLD");
  domain->step.stats   = iparam("IO_STAT");

  /* These look circular, but they're really not */

  if (domain->step.measure)
    domain->mea = measure_alloc(domain);
  if (domain->step.stats) 
    domain->stats = stat_alloc (domain);
  
  /* Allocate matrices */

  Re = dparam("Re");
  dt = dparam("DT");
  Je = iparam("TORDER");

  if (user_build) {
    Prism_log(0,"\nUser-defined build:\n");

    domain->Pressure = (*user_build)(domain->P, domain->Pbc, session, 0.);
    domain->Velocity = (*user_build)(domain->U, domain->Ubc, session, 
				     get_gamma(Je)*Re/dt);
    Prism_log (0,"\n");
  } else {
    Prism_log (0,"\n");
    Prism_log (0,"Building pressure matrices [");
    iparam_set ("MWORDS", iparam("MWORDSP"));
    domain->Pressure = buildops(domain->P, domain->Pbc, session, 0.);
    Prism_log (0,"]\n");
    
    Prism_log (0,"Building velocity matrices [");
    iparam_set ("MWORDS", iparam("MWORDSU"));
    domain->Velocity = buildops(domain->U, domain->Ubc, session, 
				get_gamma(Je) * Re/dt);
    Prism_log (0,"]\n");
  }

  if (option("prep")) {            /* Pre-processor exit */
    Prism_log(0,"\nDone.\n");
    Prism_exit();
  }

  /* Open the output files */

  sprintf (fname, "%s.fld", session);
  domain->fld_file = fopen (fname, "w");
  if (domain->his_list) {
    sprintf (fname, "%s.his", session);
    domain->his_file = fopen (fname, "w");
  }

#if DIM == 3
  Transform (domain->U, *domain->U->base, Fourier);
  Transform (domain->V, *domain->V->base, Fourier);
  Transform (domain->W, *domain->W->base, Fourier);
#endif

  return domain;
}

void Domain_free (Domain *domain)
{
  int n;
  const int Je = iparam("TORDER");

  Field_free (domain->P);
  Field_free (domain->U);
  Field_free (domain->V);
#if DIM==3
  Field_free (domain->W);
#endif

  for (n = 0; n < Je; n++) {
    Field_free (domain->Us[n]);
    Field_free (domain->Uf[n]);
    Field_free (domain->Vs[n]);
    Field_free (domain->Vf[n]);
#if DIM==3
    Field_free (domain->Ws[n]);
    Field_free (domain->Wf[n]);
#endif
  }

  BC_free (domain->Ubc);
  BC_free (domain->Vbc);
  BC_free (domain->Pbc);
#if DIM==3
  BC_free (domain->Wbc);
#endif

  if (domain->mea)
    measure_free(domain->mea);
  if (domain->stats)
    stat_free(domain->stats);

  /* free Matrices */
  /* free Green's function */
  /* free history points */
  /* free optional fields */

  Family_destroy();

  if (domain->fld_file)
    fclose (domain->fld_file);
  if (domain->his_file)
    fclose (domain->his_file);

  free (domain->name);
  free (domain);
}

#define SQR(x) ((x)*(x))


/* ---------------------------------------------------------------------- *
 * File_backup() -- Create a backup copy of a file                        *
 * ---------------------------------------------------------------------- */

extern int unlink (const char *path);
extern int link   (const char *path1, const char *path2);

static int File_backup (char *path1)
{
  char path2[FILENAME_MAX];
  int  stat;

  sprintf (path2, "%s.bak", path1);
  unlink  (path2);                    /* unlink path2 regardless    */
  if (!(stat = link(path1, path2)))   /* try to link path1 -> path2 */
    unlink (path1);                   /* unlink path1 only if the   */
  return stat;                        /* link was sucessful         */
}

void Domain_checkpt (Domain *domain)
{
  FILE *fp = domain->fld_file;
  char fname[FILENAME_MAX];

  sprintf(fname, "%s.chk", domain->name);
  ROOTONLY {
    File_backup(fname);
  }

  /* Save the checkpoint file and restore the permanent file pointer */

  Prism_log(1,"Checkpoint: ");

  GSYNC;
  domain->fld_file = fopen(fname, "w");
  Domain_save(domain);
  fclose(domain->fld_file);
  domain->fld_file = fp;

  /* flush the other output files */

  ROOTONLY {
    if (domain->his_file) fflush(domain->his_file);
    if (domain->mea)      fflush(domain->mea->fp);
  }
}

void Domain_save (Domain *domain)
{
  FieldFile *f = FieldFile_alloc();
  int i;

#if DIM==3
  Transform (domain->U, *domain->U->base, Physical);
  Transform (domain->V, *domain->V->base, Physical);
  Transform (domain->W, *domain->W->base, Physical);
  Transform (domain->P, *domain->P->base, Physical);
#endif  

  if (iparam("EQTYPE") == Rotational) {
    const int n = domain->U->nr * domain->U->ns * domain->U->nz *
      Field_count(domain->U);

    for (i = 0; i < n; i++) {
#if DIM==2
      (*domain->P->base)[i] -= 0.5 * 
	( SQR((*domain->U->base)[i]) + 
	  SQR((*domain->V->base)[i]) );
#else
      (*domain->P->base)[i] -= 0.5 * 
	( SQR((*domain->U->base)[i]) + 
	  SQR((*domain->V->base)[i]) + 
	  SQR((*domain->W->base)[i]) );
#endif
    }
  }

  Prism_log(0, "Writing field file [step=%d, time=%g]\n", 
	    FIELDFILE_STEP(f), FIELDFILE_TIME(f));

  FieldFile_setName(f,domain->name);
  FieldFile_put    (f,domain->U);
  FieldFile_put    (f,domain->V);
#if DIM==3
  FieldFile_put    (f,domain->W);
#endif
  FieldFile_put    (f,domain->P);
  FieldFile_write  (f,domain->fld_file);
  FieldFile_free   (f);

#if DIM==3
  Transform (domain->U, *domain->U->base, Fourier);
  Transform (domain->V, *domain->V->base, Fourier);
  Transform (domain->W, *domain->W->base, Fourier);
#endif  

  /* flush the other output files */

  ROOTONLY {
    if (domain->his_file) fflush(domain->his_file);
    if (domain->mea)      fflush(domain->mea->fp);
  }

  if (domain->stats) stat_write(domain->stats);
}

#undef SQR

