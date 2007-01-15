#include <dns.h>

void nonLinear (Domain*    D ,
		AuxField** Us,
		AuxField** Uf)
// ---------------------------------------------------------------------------
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u), and
// *** BODY-FORCE TERMS for the INERTIA-WAVE PROBLEM ***
//
// Here N(u) represents the nonlinear advection terms in the N--S equations
// transposed to the RHS and ff is a vector of body force per unit mass.
//
// Velocity field data areas of D and first level of Us are swapped, then
// the next stage of nonlinear forcing terms N(u) - a are computed from
// velocity fields and left in the first level of Uf.
//
// Nonlinear terms N(u) are computed in skew-symmetric form (Zang 1991)
//                 ~ ~
//           N  = -0.5 ( u . grad u + div uu )
//           ~           ~        ~       ~~
//
// i.e., in Cartesian component form
//
//           N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ).
//            i           j    i      j      i j      j
//
// in cylindrical coordinates
//
//           Nx = -0.5 {ud(u)/dx + vd(u)/dy +  d(uu)/dx + d(vu)/dy +
//                 1/y [wd(u)/dz + d(uw)/dz + vu      ]}
//           Ny = -0.5 {ud(v)/dx + vd(v)/dy +  d(uv)/dx + d(vv)/dy +
//                 1/y [wd(v)/dz + d(vw)/dz + vv - 2ww]}
//           Nz = -0.5 {ud(w)/dx + vd(w)/dy +  d(uw)/dx + d(vw)/dy +
//                 1/y [wd(w)/dz + d(ww)/dz + 3wv     ]}
//
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in reference[3].
//
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
//
// FORCING FOR INERTIA-WAVE PROBLEM: SEE BELOW.
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();  // -- Number of space dimensions.   
  const int_t NCOM = D -> nField() - 1; // -- Number of velocity components. 
  int_t       i, j;

#if defined (ALIAS)
  const int_t       nZ32   = Geometry::nZProc();
#else
  const int_t       nZ32   = Geometry::nZ32();
#endif 
  const int_t       nZ     = Geometry::nZ();
  const int_t       nZP    = Geometry::nZProc();
  const int_t       nP     = Geometry::planeSize();
  const int_t       nPP    = Geometry::nBlock();
  const int_t       nPR    = Geometry::nProc();
  const int_t       nTot   = Geometry::nTotProc();
  const int_t       nTot32 = nZ32 * nP;
  vector<real_t>    work ((2 * NCOM + 1) * nTot32);
  vector<real_t*>   u32 (NCOM), n32 (NCOM);
  vector<AuxField*> U   (NCOM), N   (NCOM);
  Field*            master = D -> u[0];
  real_t*           tmp    = &work[0] + 2 * NCOM * nTot32;

  Veclib::zero ((2 * NCOM + 1) * nTot32, &work[0], 1); // -- Catch-all cleanup.

  for (i = 0; i < NCOM; i++) {
    u32[i] = &work[0] +  i         * nTot32;
    n32[i] = &work[0] + (i + NCOM) * nTot32;

    AuxField::swapData (D -> u[i], Us[i]);

    U[i] = Us[i];
    U[i] -> transform32 (INVERSE, u32[i]);
    N[i] = Uf[i];
  }

  for (i = 0; i < NCOM; i++) {

    // -- Terms involving azimuthal derivatives and frame components.

    if (i == 0)
      Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
    if (i == 1)
      Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);

    if (NCOM == 3) {

      if (i == 1)
	Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
      if (i == 2)
	Veclib::svvtt   (nTot32,  3.0, u32[2], 1, u32[1], 1,      n32[2], 1);

      if (nZ > 2) {
	Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ, nPP, tmp, 2);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32[i], 1, n32[i], 1);

	Veclib::vmul       (nTot32, u32[i], 1, u32[2], 1, tmp, 1);
	Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ, nPP, tmp, 2);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }
    }

    if (i == 2) master -> divY (nZ32, n32[i]);

    // -- 2D non-conservative derivatives.
    
    for (j = 0; j < 2; j++) {
      Veclib::copy (nTot32, u32[i], 1, tmp, 1);
      master -> gradient (nZ32, nP, tmp, j);

      if (i <  2) master -> mulY (nZ32, tmp);

      Veclib::vvtvp (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
    }

    // -- 2D conservative derivatives.
     
    for (j = 0; j < 2; j++) {
      Veclib::vmul (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
      master -> gradient (nZ32, nP, tmp, j);

      if (i <  2) master -> mulY (nZ32, tmp);

      Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
    }

    Blas::scal (nTot32, -0.5, n32[i], 1);
  }

#if 0
  // -- Add forcing for inertia-wave problem (transient version):
  //    (NB: T here is actually OMEGA_1*t.)
  //
  // Fx: +[2*w2/w1*sin(THETA)+THETA_DDOT] *                y*sin(T-z)
  //     +[w2*THETA_DOT/(w1*w1)*cos(THETA)-THETA_DOT/w1] * y*cos(T-z)
  //     +2*{w2/w1*sin(THETA)*cos(T)-THETA_DOT/w1*sin(T)}*{v*cos(z)-w*sin(z)}
  //     +2*{w2/w1*sin(THETA)*sin(T)-THETA_DOT/w1*cos(T)}*{v*sin(z)+w*cos(z)}
  //
  // Fy: +2*[w2/w1*cos(THETA)+1]*w
  //     -2*[w2/w1*sin(THETA)*cos(T-z)+THETA_DOT/w1*sin(T-z)]*u
  //     -{[w2*THETA_DOT/(w1*w1)*cos(THETA)-THETA_DOT/w1]*cos(T-z)
  //                                       +THETA_DDOT   *sin(T-z)}*x
  //
  // Fz: -2*[w2/w1*cos(THETA)+1]*v
  //     -2*[w2/w1*sin(THETA)*sin(T-z)-THETA_DOT/w1*cos(T-z)]*u
  //     -{[w2*THETA_DOT/(w1*w1)*cos(THETA)-THETA_DOT/w1]*sin(T-z)
  //                                       -THETA_DDOT   *cos(T-z)}*x
  //
  // The model adopted for angular motion of spin axis is harmonic:
  // if t < T_RISE:
  //
  // THETA      = 0.5*      ANG_MAX*(1-cos(w3*t)), otherwise THETA = ANG_MAX
  // THETA_DOT  = 0.5*   w3*ANG_MAX*   sin(w3*t),  otherwise THETA_DOT  = 0
  // THETA_DDOT = 0.5*w3*w3*ANG_MAX*   cos(w3*t),  otherwise THETA_DDOT = 0
  //
  // where w3 = PI/T_RISE.
  //
  // NBNBNB: We still need to multiply axial and radial terms by y.

  static const real_t Tr    = Femlib::value ("T_RISE");
  static const real_t Th    = Femlib::value ("ANG_MAX");
  static const real_t w1    = Femlib::value ("OMEGA_1");
  static const real_t w2    = Femlib::value ("OMEGA_2");
  static const real_t w3    = Femlib::value ("PI/T_RISE");
  static const real_t dz    = Femlib::value ("TWOPI/BETA") / (nZ32 * nPR);
  static const real_t dt    = Femlib::value ("D_T");
  static const real_t w2ow1 = w2/w1;

  const real_t t            = Femlib::value ("t") - dt;
  const real_t T            = w1 * t;
  const int  transient    = t < Tr;

  real_t theta, thetaD, thetaDD;
  
  if (transient) {
    theta   = 0.5*Th*(1.0-cos(w3*t));
    thetaD  = 0.5*w3*Th*sin(w3*t);
    thetaDD = 0.5*w3*w3*Th*cos(w3*t);
  } else {
    theta   = Th;
    thetaD  = thetaDD = 0.0;
  }

  // -- V & w cross terms.

  const real_t cross = 2.0*w2ow1*cos(theta)+1.0;
 
  Blas::axpy (nTot32,  cross, u32[2], 1, n32[1], 1);
  Blas::axpy (nTot32, -cross, u32[1], 1, n32[2], 1);

  // -- Everything else requires z position.

  const real_t cost = cos (T);
  const real_t sint = sin (T);
  const real_t x1   =  2.0*w2ow1*sin(theta)+thetaDD;
  const real_t x2   =  w2*thetaD/(w1*w1)*cos(theta)-thetaD/w1;
  const real_t x3   =  2.0*(w2ow1*sin(theta)*cost-thetaD/w1*sint);
  const real_t x4   =  2.0*(w2ow1*sin(theta)*sint-thetaD/w1*cost);

  real_t z, costmz, sintmz;

  for (i = 0; i < nZ32; i++) {

    z      = dz * (i + nZ32 * Geometry::procID());
    costmz = cos (T - z);
    sintmz = sin (T - z);

    Veclib::fill (nP, 1.0, tmp, 1);
    master -> mulY (1, tmp);
  
    // -- Axial momentum.

    Blas::axpy (nP, x1*sintmz+x2*costmz, tmp,         1, n32[0]+i*nP, 1);
    Blas::axpy (nP, x3*cos(z)+x4*sin(z), u32[1]+i*nP, 1, n32[0]+i*nP, 1);
    Blas::axpy (nP, x4*cos(z)-x3*sin(z), u32[2]+i*nP, 1, n32[0]+i*nP, 1);

    Veclib::fill (nP, 1.0, tmp, 1);
    master -> mulX (1, tmp);

    // -- Radial momentum.

    Blas::axpy (nP, -2.0*(w2ow1*sin(theta)*costmz+thetaD/w1*sintmz),
		u32[0]+i*nP, 1, n32[1]+i*nP, 1);
    if (transient)
      Blas::axpy (nP, -x2*costmz-thetaDD*sintmz, tmp, 1, n32[1]+i*nP, 1);

    // -- Angular momentum.

    Blas::axpy (nP, -2.0*(w2ow1*sin(theta)*sintmz+thetaD/w1*costmz),
		u32[0]+i*nP, 1, n32[2]+i*nP, 1);
    if (transient)
      Blas::axpy (nP, -x2*sintmz+thetaDD*costmz, tmp, 1, n32[2]+i*nP, 1);
  }

  for (i = 0; i < NDIM; i++) {
    N[i]   -> transform32 (FORWARD, n32[i]);
    master -> smooth (N[i]);
  }

#else

  // -- Shuhei's linear case. Note that the terms are taken to the RHS
  //    of the equations (i.e. negated), and in addition the axial and
  //    radial terms are multiplied by y.

  static const real_t dz    = Femlib::value ("TWOPI/BETA") / (nZ32 * nPR);
  static const real_t dt    = Femlib::value ("D_T");
  static const real_t w2ow1 = Femlib::value ("OMEGA_2 / OMEGA_1");
  static const real_t THETA = Femlib::value ("THETA");
  const real_t t            = Femlib::value ("t") - dt;
  real_t       z, cosfac;

  // -- Axial component.
  
  for (i = 0; i < nZ32; i++) {
    z      = dz * (i + nZ32*Geometry::procID());
    cosfac = -THETA*2.0*w2ow1*cos (z + t);
    Veclib::fill   (nP, cosfac, tmp, 1);
    master -> mulY (1, tmp);
    master -> mulY (1, tmp);
    Veclib::vadd (nP, n32[0] + i*nP, 1, tmp, 1, n32[0] + i*nP, 1);
  }

  // -- Radial component.

  Veclib::smul   (nTot32,  THETA*2.0*(w2ow1+1.0), u32[2], 1, tmp, 1);
  master -> mulY (nZ32, tmp);
  Veclib::vadd   (nTot32, tmp, 1, n32[1], 1, n32[1], 1);

  // -- Azimuthal component.

  Veclib::svtvp  (nTot32, -THETA*2.0*(w2ow1+1.0), u32[1],1,n32[2],1,n32[2],1);

#endif
  
  // -- Complete everything by going back to Fourier space.
  //    Add mass-matrix smoothing to enforce continuity.

  for (i = 0; i < NCOM; i++) {
   N[i] -> transform32 (FORWARD, n32[i]);
   master -> smooth    (N[i]);
  }
}
