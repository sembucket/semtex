///////////////////////////////////////////////////////////////////////////////
// particle.C
//
// Functions for integrating positions of particles.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>

Domain* FluidParticle::D     = 0;
int     FluidParticle::NDIM  = 0;
int     FluidParticle::NEL   = 0;
int     FluidParticle::NZ    = 0;
int     FluidParticle::TORD  = 0;
real*   FluidParticle::coeff = 0;
real    FluidParticle::DT    = 0.0;
real    FluidParticle::Lz    = 0.0;


FluidParticle::FluidParticle (Domain*   d,
			      const int i,
			      Point&    p) :

                              id       (i),
                              P        (p)
// ---------------------------------------------------------------------------
// Initially particle is located at p.  Find it in the 2D mesh.  Trim to
// periodic length in 3D if required.
// ---------------------------------------------------------------------------
{
  if (!D) {
    D     = d;
    NDIM  = Geometry::nDim();
    NEL   = Geometry::nElmt();
    NZ    = Geometry::nZ();
    TORD  = (int) Femlib::value ("N_TIME");
    DT    =       Femlib::value ("D_T");
    coeff = new real [TORD + 1];
    Lz    = Femlib::value ("TWOPI / BETA");
  }

  // -- Try to locate particle, stop if can't.

  register int k;
  const int    guess = 1;

  E = 0;
  for (k = 0; k < NEL; k++) {
    r = s = 0.0;
    if (D -> Esys[k] -> locate (P.x, P.y, r, s, guess)) {
      E = D -> Esys[k];
      break;
    }
  }

  if (!E) return;

  if (NDIM == 2) {
    P.z = 0.0;

    u = new real [TORD + TORD];
    v = u + TORD;

  } else {
    if   (P.z < 0.0) P.z = Lz - fmod (fabs (P.z), Lz);
    else             P.z = fmod (P.z, Lz);
    
    u = new real [TORD + TORD + TORD];
    v = u + TORD;
    w = v + TORD;
  }
}


void FluidParticle::integrate (const int step)
// ---------------------------------------------------------------------------
// Integrate massless particle's position using predictor--corrector method.
// If particles leave 2D mesh they marked by setting E = 0.  For 3D, they
// get put back into the fundamental period of the solution if they leave
// in the z-direction.
//
// NB: The domain velocity fields must be in Fourier space prior to call (3D).
// ---------------------------------------------------------------------------
{
  if (!E) return;

  register int i;
  const int    N     = min (step, TORD);
  const int    NP    = N + 1;
  const int    guess = 1;
  real         xp, yp, zp, up, vp, wp;

  if (NDIM == 2) {		// -- 2D integration.
    
    // -- Predictor.
    
    Integration::AdamsBashforth (N, coeff);
    Blas::scal (N, DT, coeff, 1);

    u[0] = D -> u[0] -> probe (E, r, s, 0);
    v[0] = D -> u[1] -> probe (E, r, s, 0);

    xp = P.x;
    yp = P.y;
    for (i = 0; i < N; i++) {
      xp += coeff[i] * u[i];
      yp += coeff[i] * v[i];
    }

    if (!E -> locate (xp, yp, r, s)) {
      E = 0;
      for (i = 0; i < NEL; i++) {
	r = s = 0.0;
	if (D -> Esys[i] -> locate (xp, yp, r, s, guess)) {
	  E = D -> Esys[i];
	  break;
	}
      }
      if (!E) return;
    }

    // -- Corrector.

    Integration::AdamsMoulton (NP, coeff);
    Blas::scal (NP, DT, coeff, 1);

    up = D -> u[0] -> probe (E, r, s, 0);
    vp = D -> u[1] -> probe (E, r, s, 0);

    P.x += coeff[0] * up;
    P.y += coeff[0] * vp;
    for (i = 1; i < NP; i++) {
      P.x += coeff[i] * u[i - 1];
      P.y += coeff[i] * v[i - 1];
    }

    // -- Maintain multilevel storage.

    if (!E -> locate (P.x, P.y, r, s)) {
      E = 0;
      for (i = 0; i < NEL; i++) {
	r = s = 0.0;
	if (D -> Esys[i] -> locate (P.x, P.y, r, s, guess)) {
	  E = D -> Esys[i];
	  break;
	}
      }
      if (!E) return;
    }

    roll (u, TORD);
    roll (v, TORD);

  } else {			// -- 3D integration.
    
    // -- Predictor.
    
    Integration::AdamsBashforth (N, coeff);
    Blas::scal (N, DT, coeff, 1);
    
    u[0] = D -> u[0] -> probe (E, r, s, P.z);
    v[0] = D -> u[1] -> probe (E, r, s, P.z);
    w[0] = D -> u[2] -> probe (E, r, s, P.z);

    xp = P.x;
    yp = P.y;
    zp = P.z;
    for (i = 0; i < N; i++) {
      xp += coeff[i] * u[i];
      yp += coeff[i] * v[i];
      zp += coeff[i] * w[i];
    }

    if (!E -> locate (xp, yp, r, s)) {
      E = 0;
      for (i = 0; i < NEL; i++) {
	r = s = 0.0;
	if (D -> Esys[i] -> locate (xp, yp, r, s, guess)) {
	  E = D -> Esys[i];
	  break;
	}
      }
      if (!E) return;
    }

    // -- Corrector.

    Integration::AdamsMoulton (NP, coeff);
    Blas::scal (NP, DT, coeff, 1);

    up = D -> u[0] -> probe (E, r, s, zp);
    vp = D -> u[1] -> probe (E, r, s, zp);
    wp = D -> u[2] -> probe (E, r, s, zp);

    P.x += coeff[0] * up;
    P.y += coeff[0] * vp;
    P.z += coeff[0] * wp;
    for (i = 1; i < NP; i++) {
      P.x += coeff[i] * u[i - 1];
      P.y += coeff[i] * v[i - 1];
      P.z += coeff[i] * w[i - 1];
    }

    // -- Maintain multilevel storage.

    if (!E -> locate (P.x, P.y, r, s)) {
      E = 0;
      for (i = 0; i < NEL; i++) {
	r = s = 0.0;
	if (D -> Esys[i] -> locate (P.x, P.y, r, s, guess)) {
	  E = D -> Esys[i];
	  break;
	}
      }
      if (!E) return;
    }
    if   (P.z < 0.0) P.z = Lz - fmod (fabs (P.z), Lz);
    else             P.z = fmod (P.z, Lz);

    roll (u, TORD);
    roll (v, TORD);
    roll (w, TORD);
  }
}
