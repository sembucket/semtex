/*
 * floK -- Dwight Barkley's stability code for steady and time-periodic flows
 *
 * General Information
 * -------------------
 * This program computes the leading eigenspectrum for two related 
 * stability problems:
 *
 *    (1) 2D or 3D perturbations to steady flows
 *
 *    (2) 3D perturbations to time-periodic flows
 *
 * The only parameter for the perturbation is the wave number, defined in
 * the input file as BETA.  This is related to the wavelength in the normal
 * way: beta = 2*PI/lambda.  The (x,y)-structure of the eigenmode is one of 
 * the results that comes out of the stability calculation.
 *
 * Other parameters and instructions for steady or time-periodic cases are
 * defined in the relevant section below.  In either case, the eigenvalue 
 * problem is solved by subspace iteration.  In the input file, the following 
 * parameters can be set to control the search for leading eigenvalues:
 *
 *    KDIM               dimension of the subspace
 *
 *    EVTOL              eigenvalue tolerance, | u - \mu A u |
 *
 *    NVEC               number of eigenvectors to test for convergence
 *
 *    NITS               maximum number of subspace iterations
 *
 *    NWRT               number of eigenvectors to write to the output
 *                       file at the end of the iteration
 *
 * There are a few additional parameters needed for the time-periodic case,
 * i.e. the period of the base flow.  These are described in the relevant
 * section below.
 *
 * The base flow for the stability problem must be computed beforehand,
 * unless it can be defined analytically.  Use the 2D version of Prism to
 * solve for the base flow.  The stability problem is then defined on the 
 * same mesh but with different boundary conditions for the perturbations.  
 * The perturbation field must vanish on the boundaries of the domain.  
 * Derivative boundary conditions, i.e. for outflow boundaries, apply to both 
 * the base flow and the perturbations.
 *
 * 2D/3D Stability Analysis of Steady Flows
 * ----------------------------------------
 * For this problem you must supply a single base flow that represents
 * a steady solution to the Navier-Stokes equations.  The linear operator A(u)
 * is defined as the integral of a perturbation to this base flow over a 
 * certain time interval T, defined as "TIME" in the input file.  In this int-
 * erval, the perturbation will grow by a factor of exp \alpha T, where \alpha
 * is the leading eigenvalue of A(u).  Computed eigenvalues are rescaled by 
 * the integration time.
 *
 * In the main input file, define "NSLICE" = 1 to identify the base flow as
 * steady-state.  The actual base flow should be defined in a second input
 * file named "BASE.rea" which contains only "INITIAL CONDITIONS". There, 
 * you can specify the base flow either analytically or numerically.
 * 
 * Example (contents of BASE.rea)
 * ------------------------------
 * Analytically-defined base flow:
 * 3 INITIIAL CONDITIONS
 * Given
 *     u = 1-y*y
 *     v = 0
 *
 * Numerically-defined base flow:
 * 2 INITIAL CONDITIONS
 * Restart
 *     channel.fld
 *
 *
 * 3D Stability Analysis of Time-Periodic Flows
 * --------------------------------------------
 * (This part isn't complete yet)
 *
 *     NSLICE          number of field dumps in the base file
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "prism/prism.h"
#include "veclib/veclib.h"

#define Vector floK_Vector              /* Prism also has a "Vector" now */

typedef struct vctr  {                  /* ------------ VECTOR ------------- */
#if DIM==3
  Element  *U, *V, *W;                  /*       3D Velocity fields          */
#else
  Element  *U, *V;                      /*       2D Velocity fields          */
#endif
} Vector;


/* ------------------------------------------------------------------------- */

#define BIG_RESIDS  0

/*
 * Global Variables 
 * ---------------------- */

static  int          tsteps, nsteps;    /* Discrete time (step)       */
static  double       dt,      ftime;    /* Continuous time            */
static  Domain       *Omega;            /* Solution domain            */
static  int          verbose;

static  char        rcsid[] = "$Revision$"; /* from drive.c */

/*
 * Prototypes
 * ---------------------- */

int  EV_test      (int iter, int kdim, double *z_vec, double *wr, double *wi, 
		   double hkkm1, double evtol, int nvec);
void EV_post      (Vector *tvecs, Vector *evecs, int ntot, int kdim, int nwrt, 
		   double *z_vec, double *wr, double *wi, int idone);
void EV_small     (Vector *tvecs, int ntot, int kdim, double *z_vec, 
		   double *wr, double *wi, double *hkkm1);
void EV_big       (Vector *tvecs, Vector *evecs, int ntot, int kdim, int nwrt, 
		   double *z_vec, double *wr, double *wi);
void EV_sort      (double *z_vec, double *wr, double *wi, double *s_flag, 
		   int kdim);
void EV_write     (int dim, Field *vwrt[]);

int    A_op       (Vector vec);
int    a_op       (Domain *Omega, int nsteps);

void   V_scal     (int ntot, double s, Vector vec                   );
void   V_axpy     (int ntot, double s, Vector vec1,   Vector vec2   );
void   V_copy     (int ntot,           Vector vec1,   Vector vec2   );
double V_dot      (int ntot,           Vector vec1,   Vector vec2   );
double V_nrm2     (int ntot,           Vector vec                   );

void PBCs_store   (Domain *omega);
void VdgradV      (Domain *omega);

BSystem *floK_build (Field *U, Bedge *Ubc, const char *name, double lambda);

/* ------------------------------------------------------------------------- */

main(int argc, char *argv[])
{ 
  Vector       *evecs;
  Vector       *tvecs;

  double       evtol, norm, *z_vec, *wr, *wi, hkkm1;
  int          ntot, iter, nits, kdim, nvec, nwrt, i, j, idone; 

  Prism_init();
  parse_args(argc, argv);

  /* Install hooks */

#if DIM==3
  user_build = floK_build;
#endif

  Omega     = Domain_alloc(argv[argc-1]);

  verbose   = option("verbose");
  dt        = dparam("DT")     ; 
  evtol     = dparam("EVTOL")  ; 
  nsteps    = iparam("NSTEPS") ;                       /* from drive.c */
  kdim      = iparam("KDIM")   ;
  nvec      = iparam("NVEC")   ;
  nwrt      = iparam("NWRT")   ;
  nits      = iparam("NITS")   ;
  ftime     = dparam("TIME_0") ;
  tsteps    = 0;
  ntot      = FIELD_NR(Omega->U) * FIELD_NS(Omega->U) * Field_count(Omega->U);

  /* Set up storage (there are k+1 vectors)  */

  evecs = (Vector *) malloc((kdim+1)*sizeof(Vector));
  for(i=0;i<kdim+1;i++) {
    const int N = Omega->U->nr;
    const int M = Omega->U->nz;

    evecs[i].U = Field_aux (Omega->U, N, N, M, 'u');
    evecs[i].V = Field_aux (Omega->U, N, N, M, 'v');
#if DIM==3
    evecs[i].W = Field_aux (Omega->U, N, N, M, 'w');
#endif
  }

  tvecs = (Vector *) malloc((kdim+1)*sizeof(Vector));
  for(i=0;i<kdim+1;i++) {
    const int N = Omega->U->nr;
    const int M = Omega->U->nz;

    tvecs[i].U = Field_aux (Omega->U, N, N, M, 'u');
    tvecs[i].V = Field_aux (Omega->U, N, N, M, 'v');
#if DIM==3
    tvecs[i].W = Field_aux (Omega->U, N, N, M, 'w');
#endif
  }

  /* Check parameters */

  if (kdim<1)
    Prism_error ("bad parameter: KDIM must be > 1");
  if (nvec<1)
    Prism_error ("bad parameter: NVEC must be > 1");
  if (nits<kdim)
    Prism_error ("bad parameter: NITS must be >= KDIM");
  if (kdim<nvec)
    Prism_error ("bad parameter: NVEC must be <= KDIM");

  /* Parameters are OK ... echo them to the logfile */

  printf ("kdim  = %d\n", kdim);
  printf ("nvec  = %d\n", nvec);
  printf ("nwrt  = %d\n", nwrt);
  printf ("nits  = %d\n", nits);
  printf ("evtol = %g\n", evtol);
  if (verbose) printf ("ntot  = %d\n", ntot);

  /* Get starting vector and normalize */
  dcopy(ntot, *Omega->U->base, 1, *evecs[0].U->base, 1);
  dcopy(ntot, *Omega->V->base, 1, *evecs[0].V->base, 1);
#if DIM==3
  dcopy(ntot, *Omega->W->base, 1, *evecs[0].W->base, 1);
#endif
  norm = V_nrm2  (ntot, evecs[0]);
  V_scal(ntot, 1./norm, evecs[0]);
  if(verbose) printf("initial norm = %g\n\n", norm);

  /* Generate kdim additional vectors for the Krylov sequence:
   *
   *     x_0 = x
   *     x_1 = A x
   *     x_2 = A x_1 = A^2 x
   *     x_3 = A x_2 = A^3 x    ... up to x_k
   */

  for (i=1; i<=kdim; i++) {
    V_copy(ntot, evecs[i-1], evecs[i]);
    A_op  (evecs[i]);

    /* Look at what is happening as initial vectors are generated */
    if(i != kdim) {
      z_vec = dvector(0,i*i-1);
      wr    = dvector(0,i-1);
      wi    = dvector(0,i-1);

      for (j=0;j<=i;j++) 
	V_copy (ntot, evecs[j], tvecs[j]);

      EV_small (tvecs, ntot, i, z_vec, wr, wi, &hkkm1); 
      EV_test  (i, i, z_vec, wr, wi, hkkm1, evtol, i);

      free(z_vec); free(wr); free(wi);
    }
  }

  z_vec = dvector(0,kdim*kdim-1);
  wr    = dvector(0,kdim-1);
  wi    = dvector(0,kdim-1);

  /* MAIN LOOP */

  for(iter=kdim;iter<=nits;iter++) {

    if(iter!=kdim) {
      /* If not first time, shuffle vectors and act once with A_op.  *
       * Remember that everybody has to be normalized the same way   *
       * or we destroy the Krylov sequence.                          */
      norm = V_nrm2(ntot, evecs[1]);
      for(i=1; i<=kdim; i++) {
	V_scal(ntot, 1./norm,  evecs[i]);
	V_copy(ntot, evecs[i], evecs[i-1]);
      }
      A_op(evecs[kdim]);
    }

    /* Find evals and evecs in small space */
    for(i=0;i<=kdim;i++) V_copy(ntot, evecs[i], tvecs[i]);
    EV_small(tvecs, ntot, kdim, z_vec, wr, wi, &hkkm1); 


    /* Stopping test and output */
    if(idone = EV_test(iter, kdim, z_vec, wr, wi, hkkm1, evtol, nvec)) break;
  }

  /* Post Process */
  EV_post(tvecs, evecs, ntot, kdim, nwrt, z_vec, wr, wi, idone);

  return 0;   /* Normal exit */
}
/* ------------------------------------------------------------------------- */

int EV_test (int iter, int kdim, double *z_vec, double *wr, double *wi, 
	     double hkkm1, double evtol, int nvec)

{
  double re_ev, im_ev, max_resid, *resid=dvector(0,kdim-1);
  static double min_max1, min_max2;
  int i, idone;

  if(min_max1==0.) min_max1 = 1000.;
  if(min_max2==0.) min_max2 = 1000.;


  /* Sort small eigenvectors by residual */
  for(i=0;i<kdim;i++) {
    resid[i] = hkkm1*fabs(z_vec[kdim-1+i*kdim])
      /sqrt(ddot(kdim, z_vec+i*kdim, 1, z_vec+i*kdim, 1));
    if(wi[i]<0.) 
      resid[i-1] = resid[i] = hypot(resid[i-1],resid[i]);
  }    
  EV_sort(z_vec, wr, wi, resid, kdim);


  /* Stopping test */
  if(resid[nvec-1]<evtol) 
    idone = nvec;
  else if (min_max1<0.01 && resid[nvec-1]>10.*min_max1 ||
	   min_max2<0.01 && hkkm1        >10.*min_max2 )
    idone = -1;
  else
    idone = 0;

  min_max1 = MIN(min_max1,resid[nvec-1]);
  min_max2 = MIN(min_max2,hkkm1);


  /* Print useful stuff */
  printf("iter = %d:\n", iter); 
  for(i=0;i<kdim;i++){
    if(dparam("PERIOD") != 0.) {    /* print multipliers */
      re_ev = wr[i];
      im_ev = wi[i];
      printf("Mult(%d) = (%10.6g,%10.6g)             resid(%d) = %10.6g \n", 
	     i, re_ev, im_ev, i, resid[i]);
    } 
    else {                              /* print eigenvalues */
      re_ev = log(hypot(wr[i],wi[i]))/(nsteps*dt);
      im_ev = atan2(wi[i],wr[i])/(nsteps*dt);
      printf("e_val(%d) = (%10.6g,%10.6g)             resid(%d) = %10.6g \n", 
	     i, re_ev, im_ev, i, resid[i]);
    }
  }
  printf("                          |H_k,k-1| = %g \n", hkkm1);

  free(resid); 
  return(idone);
}

/*---------------------------------------------------------------------------*/

void EV_post (Vector *tvecs, Vector *evecs, int ntot, int kdim, int nwrt, 
	      double *z_vec, double *wr, double *wi, int idone)

{
  Field *VWRT[DIM];
  int i;

  if(idone==0) {
    printf("\nNot converged: Writing final vector only.\n\n");
    VWRT[0] = evecs[kdim].U;
    VWRT[1] = evecs[kdim].V;
#if DIM==3
    VWRT[2] = evecs[kdim].W;
#endif
    EV_write (DIM, VWRT);

  } else if (idone>0) {
    if (nwrt==0) {
      printf("\nConverged: ");
      printf("Writing initial vector for last KDIM iterations.\n\n");
      tsteps -= kdim*nsteps;
      ftime  -= kdim*nsteps*dt;
      VWRT[0] = evecs[0].U;
      VWRT[1] = evecs[0].V;
#if DIM==3
      VWRT[2] = evecs[0].W;
#endif
      EV_write (DIM, VWRT);

    } else {
      printf("\nConverged: Writing %d eigenvectors.\n\n", nwrt);
      EV_big(tvecs, evecs, ntot, kdim, nwrt, z_vec, wr, wi);
      for(i=0;i<nwrt;i++){
	VWRT[0] = evecs[i].U;
	VWRT[1] = evecs[i].V;
#if DIM==3
	VWRT[2] = evecs[i].W;
#endif
	EV_write (DIM, VWRT);
      }
    }
  } else {
    printf("\nMinimum residual reached: ");
    printf("Writing initial vector for last KDIM iterations.\n\n");
    tsteps -= kdim*nsteps;
    ftime  -= kdim*nsteps*dt;
    VWRT[0] = evecs[0].U;
    VWRT[1] = evecs[0].V;
#if DIM==3
    VWRT[2] = evecs[0].W;
#endif
    EV_write (DIM, VWRT);
  }
  
  fclose (Omega->fld_file);
  fclose (Omega->his_file);
}

void EV_write (int dim, Field *VWRT[])
{
  FieldFile *ff = FieldFile_alloc();
  int i;

  FieldFile_setName (ff, Omega->name);
  for (i = 0; i < dim; i++)
    FieldFile_put (ff, VWRT[i]);
  FieldFile_write (ff, Omega->fld_file);

  FieldFile_free  (ff);
}


/* ------------------------------------------------------------------------- *
 * EV_small()                                                                *
 *                                                                           *
 * Compute the approximate eigenvalues associated with a particular sub-     *
 * space.  If the iteration has actually converged to an invariant subspace, *
 * the eigenvalues will be exact.   The calculation is done in two steps:    *
 *                                                                           *
 *     (1) Given T (the Krylov sequence), find Q (an orthonormal basis).     *
 *         This is done by the Gram-Schmidt procedure which produces the     *
 *         factorization T = Q R.                                            *
 *                                                                           *
 *     (2) Compute the matrix H = Q* A Q.  If Q has converged to an          *
 *         invariant subspace, H will be a diagonal matrix containing        *
 *         the eigenvalues associated with the invariant subspace. Other-    *
 *         wise, the diagonal of H will give an approximation to the eigen-  *
 *         values and the off-diagonal terms provide an indication of the    *
 *         errors.                                                           *
 *                                                                           *
 * Input:                                                                    *
 *                                                                           *
 *         tvecs            Krylov vectors v, Av, A^2v, ...                  *
 *         ntot             dimension of the space, ntot = dim(v)            *
 *         kdim             dimension of the subspace, kdim << ntot          *
 *                                                                           *
 * Output:                                                                   *
 *                                                                           *
 *         z_vec            eigenvectors of H                                *
 *         wr               real part of the eigenvalues                     *
 *         wi               imag part of the eigenvalues                     *
 *         hkkm1            measure of how well Q has converged to an in-    *
 *                          variant subspace (eigenspace) of A               *
 *                                                                           *
 * ------------------------------------------------------------------------- */

void EV_small (Vector *tvecs, int ntot, int kdim, double *z_vec, 
	       double *wr, double *wi, double *hkkm1)

{
  int i, j, l, matz=1, ier, lwork=10*kdim;

  double* c_vec = dvector (0,(kdim+1)*(kdim+1)-1);
  double* h_vec = dvector (0, kdim   * kdim   -1);
  double* work  = dvector (0, lwork);

  dzero((kdim+1)*(kdim+1), c_vec, 1);


  /* modified Gram-Schmidt orthonormalization */

  for (i = 0; i < kdim+1; i++) {
    double gsc = V_nrm2(ntot, tvecs[i]);
    if (gsc == 0.) {
      fprintf (stderr, "basis vectors are linearly dependent\n");
      exit (-1);
    }

    c_vec[i + i*(kdim+1)] = gsc;
    V_scal (ntot, 1./gsc, tvecs[i]);
    for (j = i+1; j < kdim+1; j++) {
      gsc =  V_dot (ntot, tvecs[j], tvecs[i]);
      V_axpy (ntot, -gsc, tvecs[i], tvecs[j]);
      c_vec[i + j*(kdim+1)] = gsc;
    }
  }

  /* That completes the QR-decomposition.  We now have T = Q R     *
   * where T is the Krylov sequence, Q is an orthonormal basis for *
   * the subspace, and R is an upper-triangular matrix containing  *
   * the values from the Gram-Schmidt procedure.                   */
  
  printf ("R = \n");
  for (i = 0; i < kdim; i++) {
    for (j = 0; j < kdim; j++)
      printf ("%#10.7g ", c_vec[i+j*(kdim+1)]);
    printf ("\n");
  }

  /* Now we need to compute the matrix H = Q* A Q.  This takes    *
   * advantage of the factors already computed and stored in R,   *
   * and the fact that the basis for T comes from the Krylov      *
   * sequence <v, Av, A^2v, ...>.  Using Q and R, H can be        *
   * computed as follows:                                         *
   *                                                              *
   *         h_i,j = (q_i, A q_j)                                 *
   *                                                              *
   *               = (1/r_j,j) * ( r_i,j+1 - h_i,l r_l,j )        *
   *                                                              *
   * where the inner product is taken over l < j.                 */

  for (i = 0; i < kdim; i++) {
    for (j = 0; j < kdim; j++) {
      h_vec[i+j*kdim]  = c_vec[i+(j+1)*(kdim+1)] - 
	ddot(j, h_vec+i, kdim, c_vec+j*(kdim+1), 1);
      h_vec[i+j*kdim] /= c_vec[j+j*(kdim+1)];
    }
  }

  printf ("H = \n");
  for (i = 0; i < kdim; i++) {
    for (j = 0; j < kdim; j++)
      printf ("%#10.6f ", h_vec[i + j*kdim]);
    printf ("\n");
  }

  /* call LAPACK to compute e-vals and e-vecs of H */

  dgeev ('N', 'V', kdim, h_vec, kdim, wr, wi, NULL, 1, z_vec, 
	           kdim, work, lwork, ier);

  if (ier) printf ("ERROR IN DGEEV: ier = %d\n", ier);

  if (verbose==2) {
    for (i=0;i<kdim;i++) {
      printf ("e_vec(%d): \n", i);
      for (j=0;j<kdim;j++) {
	printf (" %d  %g \n", j+i*kdim, z_vec[j+i*kdim]);
      }
    }
  }

  /* Compute how much the of (k+1) component lies outside Q */

  *hkkm1 = fabs(c_vec[ kdim   + kdim   *(kdim+1)] /
		c_vec[(kdim-1)+(kdim-1)*(kdim+1)]);

  free(c_vec); 
  free(h_vec); 
  free(work);
  return;
}

/* ------------------------------------------------------------------------- */

void EV_big (Vector *tvecs, Vector *evecs, int ntot, int kdim, int nwrt, 
	     double *z_vec, double *wr, double *wi)

{
  double norm, wgt, resid;
  int i, j;

  /* generate big e-vectors */
  for(j=0;j<kdim;j++){
    dzero(ntot, *evecs[j].U->base, 1);
    dzero(ntot, *evecs[j].V->base, 1);
#if DIM==3
    dzero(ntot, *evecs[j].W->base, 1);
#endif
    for(i=0;i<kdim;i++){
      wgt = z_vec[i+j*kdim];
      V_axpy(ntot, wgt, tvecs[i], evecs[j]);
    }
  }

  /* normalize big e-vectors */
  for(i=0;i<kdim;i++){
    if(wi[i]==0.) {
      norm = V_nrm2(ntot, evecs[i]);
      V_scal(ntot, 1./norm, evecs[i]);
    }
    else if(wi[i] > 0.) {
      norm  = pow(V_nrm2(ntot, evecs[i]),  2.);
      norm += pow(V_nrm2(ntot, evecs[i+1]),2.);
      norm = sqrt(norm);
      V_scal(ntot, 1./norm, evecs[i]  );
      V_scal(ntot, 1./norm, evecs[i+1]);
      i++;
    }
  }

#if BIG_RESIDS
  /* Compute residuals of big vectors directly */
  for(i=0;i<nwrt;i++){
    V_copy(ntot, evecs[i], tvecs[0]);
    A_op(tvecs[0]);
    tsteps -= nsteps;
    ftime  -= nsteps*dt;
    if(wi[i]==0.) {
      V_axpy(ntot, -wr[i], evecs[i], tvecs[0]);
      resid = V_nrm2(ntot, tvecs[0]);
      printf("big resid(%d) = %g \n", i, resid);
    }
    else if(wi[i] > 0.) {
      V_axpy(ntot, -wr[i], evecs[i],   tvecs[0]);
      V_axpy(ntot,  wi[i], evecs[i+1], tvecs[0]);
      resid = V_nrm2(ntot, tvecs[0])/V_nrm2(ntot, evecs[i]);
      printf("big resid(%d) = %g \n", i, resid);
    }
    else {
      V_axpy(ntot, -wr[i], evecs[i],   tvecs[0]);
      V_axpy(ntot,  wi[i], evecs[i-1], tvecs[0]);
      resid = V_nrm2(ntot, tvecs[0])/V_nrm2(ntot, evecs[i]);
      printf("big resid(%d) = %g \n", i, resid);
    }
  }
#endif
}
/* ------------------------------------------------------------------------- */

void EV_sort (double *z_vec, double *wr, double *wi, double *s_flag, int kdim)

/* 
  Sorts the evals and evecs according to s_flag.
  Modelled after Num. Rec. straight insertion
*/
{
  double wr_temp, wi_temp, sf_temp, *z_temp=dvector(0,kdim-1);
  int i,j,n;
  
  for(j=1;j<kdim;j++){
    wr_temp = wr[j];
    wi_temp = wi[j];
    sf_temp = s_flag[j];
    for(n=0;n<kdim;n++) z_temp[n]=z_vec[n+j*kdim];
    i=j-1;
    while( i>=0 && (s_flag[i]>sf_temp) ){
      wr[i+1]     = wr[i];
      wi[i+1]     = wi[i];
      s_flag[i+1] = s_flag[i];
      for(n=0;n<kdim;n++) z_vec[n+(i+1)*kdim]=z_vec[n+i*kdim];
      i--;
    }
    wr[i+1]     = wr_temp;
    wi[i+1]     = wi_temp;
    s_flag[i+1] = sf_temp;
    for(n=0;n<kdim;n++) z_vec[n+(i+1)*kdim]=z_temp[n];
  }
  free(z_temp);
}
/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- */
/*             Most routines which depend on DIM are follow                  */
/* ------------------------------------------------------------------------- */

void V_scal(int ntot, double s, Vector vec)

{
  dscal(ntot, s, *vec.U->base, 1);
  dscal(ntot, s, *vec.V->base, 1);
#if DIM==3
  dscal(ntot, s, *vec.W->base, 1);
#endif
  return;
}

/* ------------------------------------------------------------------------- */

void V_axpy (int ntot, double s, Vector vec1, Vector vec2)

{
  daxpy(ntot, s, *vec1.U->base, 1, *vec2.U->base, 1);
  daxpy(ntot, s, *vec1.V->base, 1, *vec2.V->base, 1);
#if DIM==3
  daxpy(ntot, s, *vec1.W->base, 1, *vec2.W->base, 1);
#endif
  return;
}
/* ------------------------------------------------------------------------- */

void V_copy (int ntot, Vector vec1, Vector vec2)
{
  dcopy(ntot, *vec1.U->base, 1, *vec2.U->base, 1);
  dcopy(ntot, *vec1.V->base, 1, *vec2.V->base, 1);
#if DIM==3
  dcopy(ntot, *vec1.W->base, 1, *vec2.W->base, 1);
#endif
  return;
}

/* ------------------------------------------------------------------------- */

double V_dot (int ntot, Vector vec1, Vector vec2)

{
  double dot;
  dot  = ddot(ntot, *vec1.U->base, 1, *vec2.U->base, 1);
  dot += ddot(ntot, *vec1.V->base, 1, *vec2.V->base, 1);
#if DIM==3
  dot += ddot(ntot, *vec1.W->base, 1, *vec2.W->base, 1);
#endif
  return dot;
}

/* ------------------------------------------------------------------------- */

double V_nrm2 (int ntot, Vector vec)

{
  double norm;
  norm  = ddot(ntot, *vec.U->base, 1, *vec.U->base, 1);
  norm += ddot(ntot, *vec.V->base, 1, *vec.V->base, 1);
#if DIM==3
  norm += ddot(ntot, *vec.W->base, 1, *vec.W->base, 1);
#endif
  return sqrt(norm);
}
/* ------------------------------------------------------------------------- */

int A_op (Vector v)

{
  Omega->U = v.U;
  Omega->V = v.V;
#if DIM==3
  Omega->W = v.W;
#endif

  dparam_set ("TIME_0", 0.);

  tsteps += nsteps;
  ftime  += nsteps*dt;
  return a_op(Omega, nsteps);
}
/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- *
 * PBCs_store() -- Store curl of vorticity for the pressure b.c.'s           *
 *                                                                           *
 * This function is needed because floK computes the curl of vorticity in a  *
 * special way.  In Prism, there is a function called ComputePBCs that calc- *
 * ulates the curl of vorticity and uses it to initialize the PBCs.  Then,   *
 * when SetPBCs() is called it just adds in the nonlinear terms and does the *
 * time integration.                                                         *
 *                                                                           * 
 * In floK, the curl of vorticity is computed at the same time the nonlinear *
 * are (as it used to be done in Prism), so this function simply stores the  *
 * computed values.                                                          *
 *                                                                           *
 * NOTE: This function should only be called from DN().                      *
 * ------------------------------------------------------------------------- */

void PBCs_store(Domain *omega)
{
  Bedge   *Pbc    = omega->Pbc;
  BSystem *M      = omega->Pressure_sys;
  double   viscos = scalar("1./Re");    /* Viscous scaling */
  struct { Field *x, *y; } BC;

  BC.x = omega->Us[0];    
  BC.y = omega->Vs[0];

  /* .......... High-Order Pressure Boundary Conditions .......... */

  while (Pbc) {
    if (Pbc->type == 'F') {
      const int id      = Pbc->elmt->id;
      const int np      = Pbc->edge->np;
      const int start   = Pbc->edge->start;
      const int skip    = Pbc->edge->skip;

      double *nx     = Pbc->edge->unx,
             *ny     = Pbc->edge->uny,
             *dpdn   = Pbc->bc  .value,
             *curl_x = *BC.x[id].base + start,
             *curl_y = *BC.y[id].base + start;
      register int i;

      /* Store the boundary condition (2D or 3D w/nz = 1) */

      for (i = 0; i < np; i++)
	dpdn[i] = -viscos * 
	  (nx[i] * curl_x[i*skip] + ny[i] * curl_y[i*skip]);
    }
    Pbc = Pbc->next;
  }
  return;
}

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*  Below are routines from files other than drive.c and convective.c        */
/*  that needed modifications. All changes are explicitly marked with DB.    */
/*  On the Cray be careful because the Linker might not include these        */
/*  routines even though it says that it is. It is best to comment them      */
/*  out of the prism source files.                                           */
/*                                                                           */
/* ------------------------------------------------------------------------- */

#if DIM==3  /* Nothing below is needed for DIM=2 */

BSystem *floK_build (Field *U, Bedge *Ubc, const char *name, double lambda)
{
  double  beta;
  BSystem *B;
  int      k, mz, nz, direct, pid, nprocs;

  nz     = U->nz;
  beta   = dparam("BETA");
  direct = option("direct");
  pid    = option("procid");
  nprocs = option("nprocs");

  B = Matrix_alloc (U, Ubc, name);

  /* Do a couple checks to see which matrices are needed for this field.  *
   * The three solve levels (direct = [0,1,2]) require (0) no velocity    *
   * matrices, 1 pressure matrix; (1) no velocity matrices, nz pressure   *
   * matrices; and (2) nz velocity and pressure matrices.                 */

  if (lambda != 0. && direct < 2) {
    ROOTONLY fputs("none]\n", stdout);
    return B;
  } else if (direct < 1) nz = (pid ? 0 : 1);
  

  if (nprocs > 32 || nz > 32) ROOTONLY fputs("...please wait...", stdout);
  GSYNC;

  for (k = 0; k < nz; k += 2) {

    mz = (k + pid*nz) >> 1;                   /* Wave number in Z   */
#if 1 /* Changed by DB */
      /* Why?  Because when k = 0 we're really solving for k = 1    */
    B->constant = lambda + beta*beta;         /* Helmholtz constant */
    printf("beta = %g, lambda = %g ", beta, B->constant); 
#else
    B->constant = lambda + pow(beta*mz,2.);     /* Helmholtz constant */
#endif

    /* Try to load the matrix from a file.  If it doesn't exist then *
     * compute a new set and save it.  Handled by the frame routines */

    if (Frame_init (U, B, mz, name)) {

      Matrix_build (U, B);
      Frame_save   (U, B, mz);

      /* For small simulations echo a '*' for a banded system, a '.' for *
       * non-banded systems, and an 'f' if the matrix was succesfully    *
       * loaded from a file.  No echo for large simulations.             */

      if (nprocs > 32 || nz > 32) 
	continue;
      else
	putchar ((B->bandwidth > 0) ? '*' : '.');
    } else                                       
      putchar ('f');             /* ...Frame_init() loaded it from a file */
    fflush(stdout);
  }
  
  GSYNC; ROOTONLY fputs ("]\n", stdout);
  
  /* For 2-D simulations, the matrices always reside in memory, i.e., no  *
   * matrices are swapped from disk.  Since there is no call to LoadFrame *
   * in the main integration loop, the frames need to be loaded here once *
   * and for all.                                                         */

  if (DIM == 2) Frame_load (U, B, 0);

  return B;
}

/* ------------------------------------------------------------------------- */
     
void Field_gradz (Field* U, Field* Uz)
{
  int     ntot    = U->nr * U->ns * Field_count(U);
  double  beta    = dparam("BETA");

     /* Completely changed: gradz \equiv multiplication by -beta */

#if 0
  beta = -beta;
#endif
  dsmul(ntot, beta, *U->base, 1, *Uz->base, 1);

  return;
}

/* ------------------------------------------------------------------------- */

void Transform (Field *U, double *Fu, ACTION dir)

     /* Completely changed: Transform \equiv identity */

{
  int  ntot = U->nr * U->ns * Field_count(U),
       nz   = U->nz;

  dcopy (ntot* nz, *U->base, 1, Fu, 1);
}

void dfft1di () { fputs ("error: someone called dfft1di\n", stderr); }
void dfft1du () { fputs ("error: someone called dfft1du\n", stderr); }

/* ------------------------------------------------------------------------- */

#if 0
void Analyzer (Domain *omega, double time, int step)
{
  FILE     *fp = NULL;
  Element  *V[DIM + 1];
  char      fname[FILENAME_MAX];
  register  i;

  if ((step % 500) == 0) printf ("....%d\n", step);

  /* ..........  Field Files   ......... */
  
  if (step == 0 || step == (iparam("NSTEPS")) ) return; 

  V[0]   = omega->U;
  V[1]   = omega->V;
  V[2]   = omega->W;
  V[DIM] = omega->P;

  if (step % (iparam("IO_FLD")) == 0) {         
    FieldFile* tmp = FieldFile_tmp (omega->name, step, time, DIM + 1, V);
    fp = omega->fld_file;
    FieldFile_write (tmp, fp);
    FieldFile_free  (tmp);
  }
}
#endif

#endif

