/* 
 * SkewSymm() - Convective terms in skew-symmetric form
 * 
 * The following function computes the nonlinear advection term in the 
 * Navier--Stokes equations using skew-symmetric form.  The advection term
 * can be written as
 *
 *            n_i = -1/2 ( u_j du_i / dx_j + d(u_i u_j) / dx_j ).
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdlib.h>

#include "prism/prism.h"
#include "veclib/veclib.h"
#include "speclib/speclib.h"

void SkewSymm (Domain *domain)
{
#if DIM == 2
  SkewSymm_2D(domain);
#else
  SkewSymm_3D(domain);
#endif
}

int SkewSymm_2D (Domain *omega)
{
  Field   *Ux   =  omega->U,
          *Uy   =  omega->V,
          *Qz   =  omega->P,
          *Nx   = *omega->Uf,
          *Ny   = *omega->Vf,
          *Us   = *omega->Us,
          *Vs   = *omega->Vs;
  BSystem *M    =  omega->Velocity;
  const int ntot = Ux->nr * Ux->ns * M->elements;
  register int i;


  /* New pressure boundary conditions */

  ComputePBCs (omega);

  /* Conservative form */

  dvmul (ntot, *Ux->base, 1, *Ux->base, 1, *Us->base, 1);
  dvmul (ntot, *Uy->base, 1, *Uy->base, 1, *Vs->base, 1);

  Field_grad (Us, Nx, NULL);
  Field_grad (Vs, NULL, Ny);

  dvmul (ntot, *Ux->base, 1, *Uy->base, 1, *Us->base, 1);

  Field_grad (Us, Vs, Us);
  for (i = 0; i < ntot; i++) {
    (*Nx->base)[i] += (*Us->base)[i];
    (*Ny->base)[i] += (*Vs->base)[i];
  }

  /* Convective form */
  
  Field_grad (Ux, Us, Vs);
  for (i = 0; i < ntot; i++)
    (*Nx->base)[i] += (*Ux->base)[i] * (*Us->base)[i] 
                    + (*Uy->base)[i] * (*Vs->base)[i];

  Field_grad (Uy, Us, Vs);
  for (i = 0; i < ntot; i++)
    (*Ny->base)[i] += (*Ux->base)[i] * (*Us->base)[i] 
                    + (*Uy->base)[i] * (*Vs->base)[i];

  
  /* Rescale & Average */

  dscal (ntot, -.5, *Nx->base, 1);  Field_davg (Nx, M);
  dscal (ntot, -.5, *Ny->base, 1);  Field_davg (Ny, M);

  return 0;
}

/* The following form works for both serial and parallel codes */

int SkewSymm_3D (Domain *omega)
{
  BSystem *M = omega->Velocity;

  Field *U[3], *N[3], *tmp[3];

  const int nprocs = option("nprocs");
  const int ntot   = FIELD_NR(omega->U)*FIELD_NS(omega->U)*M->elements;
  const int nxy    = ntot + ntot % nprocs;
  const int nz     = FIELD_NZ(omega->U);
  const int npts   = nxy*nz;

  double **u = dmatrix(0, 2, 0, npts-1);
  double **n = dmatrix(0, 2, 0, npts-1);

  double *work1, *work2;

  int i, j;

  /* Set up field pointers */

  U[0] = omega->U;  N[0] = *omega->Uf;  tmp[0] = *omega->Us;
  U[1] = omega->V;  N[1] = *omega->Vf;  tmp[1] = *omega->Vs;
  U[2] = omega->W;  N[2] = *omega->Wf;  tmp[2] = *omega->Ws;

  work1 = *(N[1]->base);
  work2 = *(N[2]->base);

  /* Compute pressure boundary conditions */

  ComputePBCs (omega);

  /* Initialize */

  for (i = 0; i < DIM; i++)
    Transform (U[i], u[i], Physical);
  for (i = 0; i < DIM; i++)
    dzero (npts, n[i], 1);

  /* ----- Compute the term: n_i += u_j d(u_i) / dx_j ----- */

  for (i = 0; i < DIM; i++) {

    dcopy (npts, u[i], 1, *tmp[2]->base, 1);
    Field_grad_3D(tmp[2],  tmp[0], tmp[1], work1, work2);
    Field_gradz  (  U[i],  tmp[2]);
    Transform    (tmp[2], *tmp[2]->base, Physical);

    for (j = 0; j < DIM; j++)
      dvvtvp (npts, u[j], 1, *tmp[j]->base, 1, n[i], 1, n[i], 1);
  }

  /* ----- Compute the term: n_i += d(u_i u_j) / dx_j ----- */

  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {

      dvmul (npts, u[i], 1, u[j], 1, *tmp[2]->base, 1);

      switch (j) {
      case 0:
	Field_grad_3D (tmp[2], tmp[0], NULL, work1, work2);
	break;
      case 1:
	Field_grad_3D (tmp[2], NULL, tmp[0], work1, work2);
	break;
      case 2:
	Transform   (tmp[2], *tmp[2]->base, Fourier);
	Field_gradz (tmp[2],  tmp[0]);
	Transform   (tmp[0], *tmp[0]->base, Physical);
	break;
      }
      
      dvadd (npts, *tmp[0]->base, 1, n[i], 1, n[i], 1);
    }
  }
      
  /* Transform back to Fourier space */

  for (i = 0; i < DIM; i++) {
    dsmul (npts, -0.5, n[i], 1, *N[i]->base, 1);
    Transform     (N[i], *N[i]->base, Fourier);
    Field_davg_3D (N[i], M);
  }

  free_dmatrix (u, 0, 0);
  free_dmatrix (n, 0, 0);

  return 0;
}



