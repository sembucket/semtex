/* ------------------------------------------------------------------------- *
 * VxOmega() - Calculate the nonlinear terms in rotational form              *
 *                                                                           *
 * The following function calculates the nonlinear portion of the Navier-    *
 * Stokes equation in rotational form to minimize work.  It may be expressed *
 * as:                                                                       *
 *                                                                           *
 *     N(V) = [ v.wz - w.wy ] i + [ w.wx - u.wz ] j + [ u.wy - v.wx ] k      *
 *                                                                           *
 * where V = (u,v,w) are the three velocity components and (wx,wy,wz) are    *
 * the corresponding components of the vorticty vector.                      *
 *                                                                           *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "prism/prism.h"
#include "veclib/veclib.h"

void VxOmega (Domain *omega)
{
#if DIM == 3
# ifdef PARALLEL
     VxOmega_mp (omega);
#  else
     VxOmega_sp (omega);
# endif
#else
     VxOmega_2D (omega);
#endif
}

int VxOmega_2D (Domain *omega)
{
  Field   *Ux   =  omega->U,     *Uy   =  omega->V;
  Field   *Nx   = *omega->Uf,    *Ny   = *omega->Vf;
  Field   *Qz   =  omega->Qz;
  BSystem *M    =  omega->Pressure;
  Bedge   *Ubc  =  omega->Ubc;

  const int ntot = Ux->nr * Ux->ns * M->elements;
  int i;

  if (option("need.vorticity")) 
    Vorticity(omega);

  /* Compute the nonlinear products */

  for (i = 0; i < ntot; i++) {
    (*Nx->base)[i] =  (*Uy->base)[i] * (*Qz->base)[i];
    (*Ny->base)[i] = -(*Ux->base)[i] * (*Qz->base)[i];
  }

  ComputePBCs(omega);
  return 0;
}

#if DIM==3

/* 3D single-processor version */

int VxOmega_sp (Domain *omega)
{
  BSystem *M    = omega->Pressure;
  int     ntot  = omega->U->nr * omega->U->ns * M->elements;
  int     ntotz = omega->U->nz * ntot;
  int     nz    = omega->U->nz;
  int     npts  = ntot *  nz;
  int     nwork = ntot * (nz + 2);
  double  **u   = dmatrix(0, 2, 0, nwork);
  double  **q   = dmatrix(0, 2, 0, nwork);
  double  *work = dvector(0, nwork);
  int k;

  /* ................ Initialization ................. */

  Field  *Ux = omega->U,  *Qx = *omega->Us,  *Nx = *omega->Uf;
  Field  *Uy = omega->V,  *Qy = *omega->Vs,  *Ny = *omega->Vf;
  Field  *Uz = omega->W,  *Qz = *omega->Ws,  *Nz = *omega->Wf;

  if (option("need.vorticity")) 
    Vorticity (omega);

  Transform (Ux, u[0], Physical);    /* Transform to physical space */
  Transform (Uy, u[1], Physical);   
  Transform (Uz, u[2], Physical);  

  Transform (Qx, q[0], Physical);
  Transform (Qy, q[1], Physical);
  Transform (Qz, q[2], Physical);

  ComputePBCs (omega);

  /* .............. Nonlinear Products ............. */

  for (k = 0; k < ntotz; k++) {
    (*Nx->base)[k] = u[1][k] * q[2][k] - u[2][k] * q[1][k];
    (*Ny->base)[k] = u[2][k] * q[0][k] - u[0][k] * q[2][k];
    (*Nz->base)[k] = u[0][k] * q[1][k] - u[1][k] * q[0][k];
  }

  Transform (Nx, *Nx->base, Fourier);   /* Transform back */
  Transform (Ny, *Ny->base, Fourier);
  Transform (Nz, *Nz->base, Fourier);

#define STATISTICS
#define LES

  free_dmatrix(u, 0, 0);
  free_dmatrix(q, 0, 0);
  free(work);

  return 0;
}

#ifdef PARALLEL

/* 3D multiprocessor version */

#define TRANSFORM(u,f,dir) Transform_mp_half(u,f,dir)

int VxOmega_mp (Domain *omega)
{
  BSystem *M = omega->Velocity;

  const int pid    = comm_rank();
  const int nprocs = comm_size();
  const int nz     = FIELD_NZ(omega->U);
  const int ntot   = FIELD_NR(omega->U)*FIELD_NS(omega->U)*M->elements;
  const int nxy    = (ntot + nprocs - 1) / nprocs;
  const int ntotp  = nxy * nprocs;
  const int ntotz  = ntotp * nz;

  double  **u    = dmatrix(0, 2, 0, ntotz - 1);
  double  **q    = dmatrix(0, 2, 0, ntotz - 1);
  double  *work  = dvector      (0, ntotz - 1);

  int i, k;

  /* ................ Initialization ................. */

  Field *U[DIM], *Q[DIM], *N[DIM];

  U [0] = omega->U;  Q [0] = *omega->Us;  N [0] = *omega->Uf;
  U [1] = omega->V;  Q [1] = *omega->Vs;  N [1] = *omega->Vf;
  U [2] = omega->W;  Q [2] = *omega->Ws;  N [2] = *omega->Wf;

  if (option("need.vorticity")) Vorticity (omega);

  dzero (3*ntotz, &u[0][0], 1);
  dzero (3*ntotz, &q[0][0], 1);
  dzero (  ntotz,     work, 1);


  /* The FFT's and Exchange's are interleaved to try and overlap *
   * communications and computations as much as possible.        */

  packf     (nz, ntot, *U[0]->base, u[0]);
  Exchange  (ntotz, u[0], work, Full);
  unpackp   (nz, ntotp, nprocs, work, u[0]);
  
  for (k = 1; k < DIM; k++) {
    packf    (nz, ntot, *U[k]->base, u[k]);
    Exchange (ntotz,  u[k], work, Depart);
    TRANSFORM(U[k-1], u[k-1], Physical);
    Exchange (ntotz,  u[k], work, Arrive);
    unpackp  (nz, ntotp, nprocs, work, u[k]);
  }

  packf    (nz, ntot, *Q[0]->base, q[0]);
  Exchange (ntotz, q[0], work, Depart);
  TRANSFORM(U[2],  u[2], Physical);
  Exchange (ntotz, q[0], work, Arrive);
  unpackp  (nz, ntotp, nprocs, work, q[0]); 

  for (k = 1; k < DIM; k++) {
    packf    (nz, ntot, *Q[k]->base, q[k]);
    Exchange (ntotz,  q[k], work, Depart);
    TRANSFORM(Q[k-1], q[k-1], Physical);
    Exchange (ntotz,  q[k], work, Arrive);
    unpackp  (nz, ntotp, nprocs, work, q[k]);
  }

  TRANSFORM(Q[2], q[2], Physical);

  ComputePBCs (omega);

  /* .............. Non-linear Products ............. */

  for (i = 0; i < ntotz; i++) {
    (*N[0]->base)[i] = u[1][i] * q[2][i] - u[2][i] * q[1][i];
    (*N[1]->base)[i] = u[2][i] * q[0][i] - u[0][i] * q[2][i];
    (*N[2]->base)[i] = u[0][i] * q[1][i] - u[1][i] * q[0][i];
  }

#define STATISTICS
#define LES

  /* Transform back to Fourier space */

  TRANSFORM(N[0], work, Fourier);
  packp    (nz, ntotp, nprocs, work, *N[0]->base);
  Exchange (ntotz, *N[0]->base, u[0], Depart);

  TRANSFORM(N[1], work, Fourier);
  Exchange (ntotz, *N[0]->base, u[0], Arrive);
  unpackf  (nz, ntot, u[0], *N[0]->base);
  packp    (nz, ntotp, nprocs, work, *N[1]->base);
  Exchange (ntotz, *N[1]->base, u[1], Depart);

  TRANSFORM(N[2], work, Fourier);
  Exchange (ntotz,*N[1]->base, u[1], Arrive);
  unpackf  (nz, ntot, u[1], *N[1]->base);
  packp    (nz, ntotp, nprocs, work, *N[2]->base);
  Exchange (ntotz,*N[2]->base, u[2], Full);
  unpackf  (nz, ntot, u[2], *N[2]->base);

  free_dvector (work, 0);
  free_dmatrix (u, 0, 0);
  free_dmatrix (q, 0, 0);

  return 0;
}

#endif   /* ... PARALLEL */
#endif   /* ... DIM == 3 */



