/*
 * Operations in Fourier space
 *
 * $Id$
 * -------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>

#include "comm/comm.h"
#include "veclib/veclib.h"
#include "prism/prism.h"

/*
 * Compute the derivative in the z-direction
 */

void Field_gradz (Field *U, Field *Uz)
{
  const int ntot  = FIELD_NR(U)*FIELD_NS(U)*FIELD_NELMT(U);
  const int ntot2 = ntot * 2;
  const int nz    = FIELD_NZ(U);
  const int pid   = option("procid");

  double  beta    = dparam("BETA");
  double *Ur      = *U ->base;
  double *Ui      = *U ->base + ntot;
  double *dUdz_r  = *Uz->base;
  double *dUdz_i  = *Uz->base + ntot;

  int i, k = 0;

  /* zero the first two planes */

  if (pid == 0) {
    dzero(ntot2, dUdz_r, 1);
    Ur += ntot2; dUdz_r += ntot2;
    Ui += ntot2; dUdz_i += ntot2;
    k  += 2;
  }

  /* compute: dU/dz = i beta k ( Ur + i Ui) = beta k ( -Ui + i Ur ) */

  while (k < nz) {
    const double beta_k = beta * ((k + pid*nz)>>1);
    
    for (i = 0; i < ntot; i++) {     /* must work if U and Uz overlap */
      double tmp =  beta_k * Ur[i];
      dUdz_r[i]  = -beta_k * Ui[i];
      dUdz_i[i]  =  tmp;
    }
    
    k      += 2;
    Ur     += ntot2;
    dUdz_r += ntot2;
    Ui     += ntot2;
    dUdz_i += ntot2;
  }
}

/* Same as Field_gradz() but for a single element */

void Element_gradz (Element *U, Element *Uz, int ntot)
{
  const int nrns  = FIELD_NR(U)*FIELD_NS(U);
  const int nz    = FIELD_NZ(U);
  const int ntot2 = ntot * 2;
  const int pid   = option("procid");
  
  double  beta    = dparam("BETA");
  double *Ur      = *U ->base;
  double *Ui      = *U ->base + ntot;
  double *dUdz_r  = *Uz->base;
  double *dUdz_i  = *Uz->base + ntot;

  int i, k = 0;

  /* zero the first two planes */

  if (pid == 0) {
    dzero (nrns, dUdz_r, 1);
    dzero (nrns, dUdz_r +  ntot, 1);
    Ur += ntot2; dUdz_r += ntot2;
    Ui += ntot2; dUdz_i += ntot2;
    k  += 2;
  }

  /* compute: dU/dz = i beta k ( Ur + i Ui) = beta k ( -Ui + i Ur ) */

  while (k < nz) {
    const double beta_k = beta * ((k + pid*nz)>>1);

    for (i = 0; i < nrns; i++) {     /* must work if U and Uz overlap */
      double tmp =  beta_k * Ur[i];
      dUdz_r[i]  = -beta_k * Ui[i];
      dUdz_i[i]  =  tmp;
    }

    k      += 2;
    Ur     += ntot2;
    dUdz_r += ntot2;
    Ui     += ntot2;
    dUdz_i += ntot2; 
  }
}

/* ------------------------------------------------------------------------- *
 * Field_grad_3D() - Scalar gradient [3-D version]                           *
 *                                                                           *
 * This is a special 3-D version of the routine Field_grad().  It vectorizes *
 * over all frames in a mesh and so is more than just a simple loop over the *
 * frames.  This function will work correctly if the input field represents  *
 * data in either physical or Fourier space.                                 *
 *                                                                           *
 * A single scalar field may be used for both input and output if only one   *
 * component is needed.                                                      *
 * ------------------------------------------------------------------------- */

void Field_grad_3D (Field *U, Field *Ux, Field *Uy, double *Ur, double *Us)
{
  const int nr   = FIELD_NR(U);
  const int ns   = FIELD_NS(U);
  const int nz   = FIELD_NZ(U);
  const int nrns = nr * ns;
  const int nel  = FIELD_NELMT(U);
  const int ntot = nrns * nel;
  const int pid  = option("procid");

  double  **dr, **ds, *fld;
  int i, k, m;

  getops (nr, 0, 0, &dr, 0);
  getops (ns, 0, 0, &ds, 0);

  /* Compute (r,s) partials */

  dgemm ('T','N', nr, ns*nel*nz, nr, 1., *dr, nr,fld = *U->base,nr, 0.,Ur,nr);
  for (k = 0; k < nel*nz; k++, fld += nrns)
    dgemm ('N','N', nr,ns,ns, 1.,fld,nr,*ds,ns,0.,Us+k*nrns,nr);

  /* Compute -X- derivatives */

  if (Ux)
    for (k = 0; k < nel; k++)
      if (Ux[k].sx)
	for (m = 0; m < nz; m++)
	  for (i = 0; i < nrns; i++)
	    (*Ux[k].base + m*ntot)[i] = 
	    (Ur + k*nrns + m*ntot)[i] * (*Ux[k].rx)[i] +
	    (Us + k*nrns + m*ntot)[i] * (*Ux[k].sx)[i];
      else
	for (m = 0; m < nz; m++)
	  for (i = 0; i < nrns; i++)
	    (*Ux[k].base + m*ntot)[i] = 
	    (Ur + k*nrns + m*ntot)[i] * (*Ux[k].rx)[i];

  
  /* Compute -Y- derivatives */
  
  if (Uy)
    for (k = 0; k < nel; k++)
      if (Uy[k].ry)
	for (m = 0; m < nz; m++)
	  for (i = 0; i < nrns; i++)
	    (*Uy[k].base + m*ntot)[i] = 
	    (Ur + k*nrns + m*ntot)[i] * (*Uy[k].ry)[i] +
	    (Us + k*nrns + m*ntot)[i] * (*Uy[k].sy)[i];
      else
	for (m = 0; m < nz; m++)
	  for (i = 0; i < nrns; i++)
	    (*Uy[k].base + m*ntot)[i] = 
	    (Us + k*nrns + m*ntot)[i] * (*Uy[k].sy)[i];
}

/* Same as Field_grad_3d(), but for a single element and without the   *
 * extra workspace.  The argument "ntot" is the number of points in an *
 * (xy)-plane, or the jump in memory between z-planes.                 */

void Element_grad_3D (Element *U, Element *Ux, Element *Uy, int ntot)
{
  const int nr   = FIELD_NR(U);
  const int ns   = FIELD_NS(U);
  const int nz   = FIELD_NZ(U);
  const int nrns = nr * ns;
  const int pid  = option("procid");

  double **dr, **ds;
  int i, m;

  double Ur[_MAX_NORDER*_MAX_NORDER];
  double Us[_MAX_NORDER*_MAX_NORDER];
  
  getops (nr, 0, 0, &dr, 0);
  getops (ns, 0, 0, &ds, 0);
  
  for (m = 0; m < nz; m++) {
    double *u  =      *U ->base + m*ntot;
    double *ux = Ux ? *Ux->base + m*ntot : NULL;
    double *uy = Uy ? *Uy->base + m*ntot : NULL;

    /* Compute (r,s) partials */

    dgemm ('T','N', nr,ns,nr, 1., *dr, nr,  u, nr, 0., Ur, nr);
    dgemm ('N','N', nr,ns,ns, 1.,  u,  nr, *ds,ns, 0., Us, nr);

    /* Compute -X- derivatives */

    if (Ux)
      if (Ux->sx)
	for (i = 0; i < nrns; i++)
	  ux[i] = Ur[i] * (*Ux->rx)[i] + Us[i] * (*Ux->sx)[i];
      else
	for (i = 0; i < nrns; i++)
	  ux[i] = Ur[i] * (*Ux->rx)[i];
  
  
  /* Compute -Y- derivatives */
  
    if (Uy)
      if (Uy->ry)
	for (i = 0; i < nrns; i++)
	  uy[i] = Ur[i] * (*Uy->ry)[i] + Us[i] * (*Uy->sy)[i];
      else
	for (i = 0; i < nrns; i++)
	  uy[i] = Us[i] * (*Uy->sy)[i];
  }
}

/* 
 * Direct Stiffness Summation - Mass matrix averaging [3D Version]
 *
 * If Prism is compiled with the nonconforming option, then Field_davg_3D
 * turns into a call to qssum3d.  This only needs to be compiled of the
 * macro version isn't defined, i.e. for a strictly conforming code.
 */

#ifndef Field_davg_3D

void Field_davg_3D (Field *U, BSystem *M)
{
  const int pid  = option("procid");
  const int NB   = M->bpts;
  const int nel  = M->elements;
  const int nz   = FIELD_NZ(U);
  const int nb   =(FIELD_NR(U) + FIELD_NS(U) - 2) << 1;
  const int nrns = FIELD_NR(U) * FIELD_NS(U);

  int    **bmap = M->bmap;
  double *field = *U->base;
  int i, k, m;

  tempVector (valu, NB);

  for (m = 0; m < nz; m++) {
    dzero(NB, valu, 1);

    for (k = 0; k < nel; k++, field += nrns)  
      for (i = 0; i < nb; i++) 
	valu [bmap[k][i]] += (*U[k].mass)[U[k].emap[i]] * field[U[k].emap[i]];
  
    dvmul(NB, valu, 1, M->massinv, 1, valu, 1);

    field -= nrns*nel;
    for (k = 0; k < nel; k++, field += nrns) 
      for (i = 0; i < nb; i++)
	field [U[k].emap[i]] = valu [bmap[k][i]];
  }

  freeVector (valu);
}

#endif

/* ------------------------------------------------------------------------ *
 * Transform() - Fourier Transform/Inverse of a Vector Field                *
 *                                                                          *
 * This function moves a vector field from Fourier <-> Physcial space, dep- *
 * ending on the argument "dir".  The transform is assumed to be real half- *
 * symmetric with the highest mode being identically zero.  The input vec-  *
 * tor field "U" is always the source, and the vector "Fu" is always the    *
 * transformation.                                                          *
 * ------------------------------------------------------------------------ */

static int irev  = 1;
static int ireal = 1;
static int isign = 0;
static int ierr  = 0;

static double *factor = NULL;
static int    *bitrev = NULL;

void Transform_init (int nz)
{
  const int nmodes = nz >> 1;

  if (factor) {
    free (factor);
    free (bitrev);
  }

  isign  = 1;
  factor = (double*) malloc(6*nmodes*sizeof(double));
  bitrev = (int*)    malloc(6*nmodes*sizeof(int));

  fftdf (NULL, nmodes, 1, 1, 1, 0, factor, irev, bitrev, ierr, ireal);
}

void Transform (Field *U, double *Fu, ACTION dir)
{
  /* Check for initialization */

  if (isign == 0) Transform_init(iparam("NZ"));

#ifdef PARALLEL
  Transform_mp (U, Fu, dir);
#else
  Transform_sp (U, Fu, dir);
#endif  
}

/* Single-processor version */

int Transform_sp (Field *U, double *Fu, ACTION dir)
{
  int nz     = U->nz;
  int ntot   = U->nr * U->ns * Field_count(U);
  int n      = nz >> 1;
  int nskip  = 1;
  int m      = 1;
  int mskip  = 1;

  double *data = U->base[0];
  double *output = Fu;

  int i;
  tempVector (tmp, nz+2);
  
  switch (dir) {
  case Physical:
    isign   = 1;
    tmp[nz] = tmp[nz+1] = 0.;
    dzero (ntot, data + ntot, 1);
    for (i = 0; i < ntot; i++, data++, output++) {
      dcopy(nz, data, ntot, tmp, 1);
      fftdf(tmp, n, nskip, m, mskip, isign, factor, irev, bitrev, ierr, ireal);
      dcopy(nz, tmp, 1, output, ntot);
    }
    break;

  case Fourier:
    isign  = -1;
    for (i = 0; i < ntot; i++, data++, output++) {
      dcopy(nz, data, ntot, tmp, 1);
      fftdf(tmp, n, nskip, m, mskip, isign, factor, irev, bitrev, ierr, ireal);
      dsmul(nz, 1./(4.*n), tmp, 1, output, ntot);
    }
    break;

  default:
    Prism_error("Prism: invalid direction in Transform\n");
    break;
  }

  freeVector (tmp);
  return 0;
}

/* ------------------------------------------------------------------------ *
 * Transform32() - special Transform for dealiasing                         *
 *                                                                          *
 * This is a special transform for de-aliasing the non-linear terms.  It    *
 * does not allocate workspace so the input or output field (on the Cray)   *
 * must be large enough to perform the transforms.  This function operates  *
 * differently than Transform()...for a forward transform (dir = Physical)  *
 * the INPUT is from U[] and the OUTPUT is to Fu[].  For an inverse trans-  *
 * form this is reversed (Fu -> U).                                         *
 *                                                                          *
 * size(Fu) = double[ ntot * (nz * 3/2 + 2) ]                               *
 *                                                                          *
 *                                                                          *
 * NOTE: This version is just a "stub"                                      *
 * ------------------------------------------------------------------------ */

void Transform32 (Field *U, double *Fu, ACTION dir)
{
  int nz     = U->nz;
  int nz32   = nz * 3 / 2;
  int ntot   = U->nr * U->ns * Field_count(U);
  int n      = nz >> 1;
  int nskip  = 1;
  int m      = ntot;
  int mskip  = (nz32 + 2) >> 1;

  double *data   = U->base[0];
  double *output = Fu;

  int i, ierr;

  if (isign == 0) 
    Prism_error("Prism: Transform must be called before Transform32\n");
  
  switch (dir) {

  case Physical:
    isign  = 1;
    dzero (ntot, data + ntot, 1);
    for (i = 0; i < ntot; i++, data++, output += (nz32 + 2)) {
      dcopy(nz, data, ntot, output, 1);
      dzero(nz32 + 2 - nz,  output + nz, 1);
    }

    fftdf(Fu, n, nskip, m, mskip, isign, factor, irev, bitrev, ierr, ireal);
    break;

  case Fourier:
    isign = -1;
    dscal(ntot * (nz32 + 2), 1./(4.*n), output, 1);
    fftdf(output,n, nskip, m, mskip, isign, factor, irev, bitrev, ierr, ireal);
    for (i = 0; i < ntot; i++, data++, output += (nz32 + 2)) 
      dcopy(nz, output, 1, data, ntot);
    break;

  default:
    Prism_error("Prism: invalid direction in Transform32\n");
    break;
  }
}

#ifdef PARALLEL

/* ------------------------------------------------------------------------ *
 * Packing and Unpacking of data for the transpose                          *
 *                                                                          *
 * The following routines rearrange data for the different calculation      *
 * phases in the parallel computation of the nonlinear terms.               *
 *                                                                          *
 *     packf   --   pack Fourier modes                                      *
 *   unpackf   --   unpack Fourier modes                                    *
 *     packp   --   pack Physical planes                                    *
 *   unpackp   --   unpack Physical planes                                  *
 * ------------------------------------------------------------------------ */

void packf (int nz, int nxy, double *a, double *ap) {
  int i;
  for (i = 0; i < nz; i += 2, a += nxy*2, ap += 2)
    zvcmplx (nxy, a, 1, a + nxy, 1, (zcomplex*) ap, nz/2);
}

void unpackf (int nz, int nxy, double *ap, double *a) {
  int i;
  for (i = 0; i < nz; i += 2, a += nxy*2, ap += 2) {
    zvreal (nxy, (zcomplex*) ap, nz/2, a, 1);
    zvimag (nxy, (zcomplex*) ap, nz/2, a + nxy, 1);
  }
}

void packp (int nz, int np, int nprocs, double *a, double *ap) {
  const int NZ    = nz * nprocs;
  const int nxy   = np/nprocs;
  const int block = nxy*nz;
  int i, j;
  for (i = 0; i < nxy; i++)
    for (j = 0; j < nprocs; j++)
      dcopy (nz, a + i*NZ + j*nz, 1, ap + i*nz + j*block, 1);
}

void unpackp (int nz, int np, int nprocs, double *ap, double *a) {
  const int NZ    = nz * nprocs;
  const int nxy   = np/nprocs;
  const int block = nxy*nz;
  int i, j;

  for (j = 0; j < nprocs; j++)
    for (i = 0; i < nxy; i++, ap += nz)
      dcopy (nz, ap, 1, a +i*NZ + j*nz, 1);
}

/* ------------------------------------------------------------------------ *
 * Exchange() - Perform a Global Exchange                                   *
 *                                                                          *
 * This function exchanges data across the processor network.  The inputs   *
 * are the number of points "npts", the input data vector "a" and the       *
 * result vector "at" which MUST NOT be the same as "a".                    *
 *                                                                          *
 * Uses the optimal exchange algorithm for a hypercube.                     *
 * ------------------------------------------------------------------------ */

#define MSGTYPE(i) (1000+i)
#define MAXPROCS    128

void Exchange (int npts, double *a, double *at, ACTION mode)
{
  const int nprocs     = comm_size();
  const int procid     = comm_rank();
  const int block_size = npts / nprocs;
  const int msg_size   = block_size * sizeof(double);

  int i;

  static int recvtag[MAXPROCS];   /* asynchronous communication tags */
  static int sendtag[MAXPROCS];

  /* Quick error checking */

  if (a == at)
    Prism_error ("Exchange: A and trans(A) must be allocated separately");
  if (npts % nprocs) 
    Prism_error ("Exchange: npts should be a multiple of nprocs");
  if (nprocs > MAXPROCS)
    Prism_error ("Exchange: too many processors!");

  /* Copy the local block directly */

  dcopy (block_size, a + procid*block_size, 1, at + procid*block_size, 1);

  /* Exchange blocks with the other processors */

  switch (mode) {

    /* A "Full" exchange performs both send's and recv's and doesn't return *
     * until all of the communications are complete.                        */

  case Full:
    for (i = 1; i < nprocs; i++) {
      const int dest = i ^ procid;
      comm_xchg (MSGTYPE(i), a  + dest*block_size, at + dest*block_size, 
		 msg_size, dest);

      Prism_log (9, "\tExchange: stage %d complete\n", dest);
    }
    break;


    /* A "Depart" exchange schedules asynchronous send's and recv's but    *
     * does not wait for any communications to finish.  Don't mess with    *
     * the data until calling for an Arrive!                               */
       
  case Depart:

    /* Schedule asynchronous receives */

    for (i = 1; i < nprocs; i++) {
      const int dest = i ^ procid;
      recvtag[i] = comm_recvx (MSGTYPE(i), at+dest*block_size, msg_size);
    }
    comm_sync();  /* guarantee no deadlock on unbuffered sends */
    Prism_log (9,"\tExchange: recv's scheduled\n");

    /* Schedule asynchronous sends */

    for (i = 1; i < nprocs; i++) {
      const int dest = i ^ procid;
      sendtag[i] = comm_sendx (MSGTYPE(i), a+dest*block_size, msg_size, dest);
    }
    Prism_log (9,"\tExchange: send's scheduled\n");

    break;


    /* An "Arrive" is the other end of a "Depart" -- it just waits for    *
     * the asynchronous communications to finish.                         */

  case Arrive:

    for (i = 1; i < nprocs; i++)
      comm_wait (sendtag[i]);
    Prism_log (9,"\tExchange: send's complete\n");

    for (i = 1; i < nprocs; i++)
      comm_wait (recvtag[i]);
    Prism_log (9,"\tExchange: recv's complete\n");
    
    break;
    
  default:
    Prism_error ("Exchange: invalid transpose mode");
  }

  Prism_log (9,"\tExchange: done\n");
}

#undef MSGTYPE
#undef MAXPROCS

/* Multiprocessor Transform routines */

int Transform_mp (Field *U, double *Fu, FLAG dir)
{
  int nz = iparam("NZ");
  int m  = 1;

  const int pid    = comm_rank();
  const int nprocs = comm_size();
  const int ntot   = FIELD_NR(U)*FIELD_NS(U)*FIELD_NELMT(U);
  const int n      = nz >> 1;
  const int nskip  = 1;
  const int mskip  = n;

  double *input  = *U->base;
  double *output = Fu;
  double *buf    = NULL;
  double *tmp    = NULL;
  double  scal   = 0.25 / n;

  int nxy, ntotp, ntotz, i;

  nxy    = (ntot + nprocs - 1) / nprocs;
  m      =  ntot - pid * nxy;
  m      =  CLAMP (m, 0, nxy);
  ntotp  =  nxy * nprocs;
  ntotz  =  nxy * nz;
  nz     =  FIELD_NZ(U);
  buf    =  dvector (0, ntotz);
  tmp    =  dvector (0, 2*n+1);

  dzero (ntotz, buf, 1);

  Prism_log (9,"Transform_mp: packing\n");

  packf   (nz, ntot, input, buf);
  Exchange(ntotz, buf, Fu, Full);
  unpackp (nz, ntotp, nprocs, Fu, buf);

  Prism_log (9,"Transform_mp: transpose complete\n");

  switch (dir) {
  case Physical:
    isign =  1;
    tmp[2*n] = tmp[2*n+1] = 0.;
    dzero (m, buf + 1, nz*nprocs);
    for (i = 0; i < m; i++) {
      dcopy (2*n, buf + 2*i*mskip, 1, tmp, 1);
      fftdf (tmp, n, 1, 1, 1, isign, factor, irev, bitrev, ierr, ireal);
      dcopy (2*n, tmp, 1, buf + 2*i*mskip, 1);
    }
    break;

  case Fourier:
    isign = -1;
    for (i = 0; i < m; i++) {
      dcopy (2*n, buf + 2*i*mskip, 1, tmp, 1);
      fftdf (tmp, n, 1, 1, 1, isign, factor, irev, bitrev, ierr, ireal);
      dsmul (2*n, 1./(4.*n), tmp, 1, buf + 2*i*mskip, 1);
    }
    break;

  default:
    Prism_error ("Transform_mp: invalid transform direction");
    break;
  }

  Prism_log (9,"Transform_mp: fft complete\n");
  Prism_log (9,"Transform_mp: unpacking\n");

  packp   (nz, ntotp, nprocs, buf, Fu);
  Exchange(ntotz, Fu, buf, Full);
  unpackf (nz, ntot, buf, output);

  Prism_log (9,"Transform_mp: done\n");

  free (buf);
  free (tmp);
  return 0;
}

/* This is Transform_mp w/out transposing the data */

int Transform_mp_half (Field *U, double *Fu, ACTION dir)
{
  const int pid    = comm_rank();
  const int nprocs = comm_size();
  const int nz     = iparam("NZ");
  const int n      = nz >> 1;
  const int nskip  = 1;
  const int ntot   = U->nr * U->ns * Field_count(U);
  const int nxy    = (ntot + nprocs - 1) / nprocs;
  const int mskip  = n;

  int m = ntot - pid*nxy;
  double scal = 0.25/n;
  int i;

  m = CLAMP(m, 0, nxy);

  switch (dir) {
  case Fourier:
    isign = -1;
    dcopy (m*nz, *U->base, 1, Fu, 1);
    break;
  
  case Physical:
    isign =  1;
    dzero (m, Fu + 1, nz);
    break;

  default:
    Prism_error("Transform_mp_half: invalid transform direction");
    break;
  }
      
  fftdf (Fu, n, nskip, m, mskip, isign, factor,
	 irev, bitrev, ierr, ireal);
  
  if (dir == Fourier) 
    dscal (m*nz, scal, Fu, 1);
  return 0;
}

#endif
