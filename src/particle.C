///////////////////////////////////////////////////////////////////////////////
// particle.C: functions for integrating positions of particles.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>

Domain*  FluidParticle::D       = 0;
integer  FluidParticle::NDIM    = 0;
integer  FluidParticle::NEL     = 0;
integer  FluidParticle::NZ      = 0;
integer  FluidParticle::TORD    = 0;
integer  FluidParticle::ID_MAX  = 0;
real*    FluidParticle::P_coeff = 0;
real*    FluidParticle::C_coeff = 0;
real     FluidParticle::DT      = 0.0;
real     FluidParticle::Lz      = 0.0;


FluidParticle::FluidParticle (Domain*       d,
			      const integer i,
			      Point&        p) :

                              id           (i),
                              P            (p)
// ---------------------------------------------------------------------------
// Initially particle is located at p.  Find it in the 2D mesh.  Trim to
// periodic length in 3D if required.
// ---------------------------------------------------------------------------
{
  if (!D) {
    D       = d;
    NDIM    = Geometry::nDim();
    NEL     = Geometry::nElmt();
    NZ      = Geometry::nZ();
    TORD    = (integer) Femlib::value ("N_TIME");
    DT      =           Femlib::value ("D_T");
    P_coeff = new real [(size_t) TORD];
    C_coeff = new real [(size_t) (TORD + 1)];
    Lz      = Femlib::value ("TWOPI / BETA");
  }

  // -- Try to locate particle, stop if can't.

  register integer k;
  const integer    guess = 1;

  E = 0;
  for (k = 0; k < NEL; k++) {
    r = s = 0.0;
    if (D -> elmt[k] -> locate (P.x, P.y, r, s, guess)) {
      E = D -> elmt[k];
      break;
    }
  }

  if (!E) return;
  if (id > ID_MAX) ID_MAX = id;

  if (NDIM == 2) {
    P.z = 0.0;

    u = new real [(size_t) (TORD + TORD)];
    v = u + TORD;

  } else {
    if   (P.z < 0.0) P.z = Lz - fmod (fabs (P.z), Lz);
    else             P.z = fmod (P.z, Lz);
    
    u = new real [(size_t) (TORD + TORD + TORD)];
    v = u + TORD;
    w = v + TORD;
  }
}


void FluidParticle::integrate (const integer step)
// ---------------------------------------------------------------------------
// Integrate massless particle's position using predictor--corrector scheme.
// If particles leave 2D mesh they marked by setting E = 0.  For 3D, they
// get put back into the fundamental period of the solution if they leave
// in the z-direction.
//
// NB: The domain velocity fields must be in Fourier space prior to call (3D).
// ---------------------------------------------------------------------------
{
  if (!E) return;

  register integer i;
  const integer    N     = min (step, TORD);
  const integer    NP    = N + 1;
  const integer    guess = 1;
  real             xp, yp, zp, up, vp, wp;

  if (N <= TORD) {
    Integration::AdamsBashforth (N,      P_coeff   );
    Integration::AdamsMoulton   (NP,     C_coeff   );
    Blas::scal                  (N,  DT, P_coeff, 1);
    Blas::scal                  (NP, DT, C_coeff, 1);
  }

  if (NDIM == 2) {		// -- 2D integration.
    
    // -- Predictor.

    u[0] = D -> u[0] -> probe (E, r, s, (integer) 0);
    v[0] = D -> u[1] -> probe (E, r, s, (integer) 0);

    xp = P.x;
    yp = P.y;
    for (i = 0; i < N; i++) {
      xp += P_coeff[i] * u[i];
      yp += P_coeff[i] * v[i];
    }

    if (!E -> locate (xp, yp, r, s)) {
      E = 0;
      for (i = 0; i < NEL; i++) {
	r = s = 0.0;
	if (D -> elmt[i] -> locate (xp, yp, r, s, guess)) {
	  E = D -> elmt[i];
	  break;
	}
      }
      if (!E) {
#if defined(DEBUG)
	char     str[StrMax];
	sprintf (str, "Particle %1d at (%f, %f, %f) left mesh",
		 id, P.x, P.y, P.z);
	message (routine, str, WARNING);
#endif
	return;
      }
    }

    // -- Corrector.

    up = D -> u[0] -> probe (E, r, s, (integer) 0);
    vp = D -> u[1] -> probe (E, r, s, (integer) 0);

    P.x += C_coeff[0] * up;
    P.y += C_coeff[0] * vp;
    for (i = 1; i < NP; i++) {
      P.x += C_coeff[i] * u[i - 1];
      P.y += C_coeff[i] * v[i - 1];
    }

    if (!E -> locate (P.x, P.y, r, s)) {
      E = 0;
      for (i = 0; i < NEL; i++) {
	r = s = 0.0;
	if (D -> elmt[i] -> locate (P.x, P.y, r, s, guess)) {
	  E = D -> elmt[i];
	  break;
	}
      }
      if (!E) {
#if defined(DEBUG)
	char     str[StrMax];
	sprintf (str, "Particle %1d at (%f, %f, %f) left mesh",
		 id, P.x, P.y, P.z);
	message (routine, str, WARNING);
#endif
	return;
      }
    }

    // -- Maintain multilevel storage.

    rollv (u, TORD);
    rollv (v, TORD);

  } else {			// -- 3D integration.
    
    // -- Predictor.

    u[0] = D -> u[0] -> probe (E, r, s, P.z);
    v[0] = D -> u[1] -> probe (E, r, s, P.z);
    w[0] = D -> u[2] -> probe (E, r, s, P.z);

    xp = P.x;
    yp = P.y;
    zp = P.z;
    for (i = 0; i < N; i++) {
      xp += P_coeff[i] * u[i];
      yp += P_coeff[i] * v[i];
      zp += P_coeff[i] * w[i];
    }

    if (!E -> locate (xp, yp, r, s)) {
      E = 0;
      for (i = 0; i < NEL; i++) {
	r = s = 0.0;
	if (D -> elmt[i] -> locate (xp, yp, r, s, guess)) {
	  E = D -> elmt[i];
	  break;
	}
      }
      if (!E) return;
    }

    // -- Corrector.

    up = D -> u[0] -> probe (E, r, s, zp);
    vp = D -> u[1] -> probe (E, r, s, zp);
    wp = D -> u[2] -> probe (E, r, s, zp);

    P.x += C_coeff[0] * up;
    P.y += C_coeff[0] * vp;
    P.z += C_coeff[0] * wp;
    for (i = 1; i < NP; i++) {
      P.x += C_coeff[i] * u[i - 1];
      P.y += C_coeff[i] * v[i - 1];
      P.z += C_coeff[i] * w[i - 1];
    }

    if (!E -> locate (P.x, P.y, r, s)) {
      E = 0;
      for (i = 0; i < NEL; i++) {
	r = s = 0.0;
	if (D -> elmt[i] -> locate (P.x, P.y, r, s, guess)) {
	  E = D -> elmt[i];
	  break;
	}
      }
      if (!E) return;
    }
    if   (P.z < 0.0) P.z = Lz - fmod (fabs (P.z), Lz);
    else             P.z = fmod (P.z, Lz);

    // -- Maintain multilevel storage.

    rollv (u, TORD);
    rollv (v, TORD);
    rollv (w, TORD);
  }
}
