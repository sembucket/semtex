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
// ---------------------------------------------------------------------------
// Initially particle is located at p.  Find it in the 2D mesh.  Trim to
// periodic length in 3D if required.
// ---------------------------------------------------------------------------
  _id   (i),
  _step (0),
  _p    (p)
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

  _E = 0;
  for (k = 0; k < NEL; k++) {
    _r = _s = 0.0;
    if (D -> elmt[k] -> locate (_p.x, _p.y, _r, _s, guess)) {
      _E = D -> elmt[k];
      break;
    }
  }

  if (!_E) return;
  if (_id > ID_MAX) ID_MAX = _id;

  if (NDIM == 2) {
    _p.z = 0.0;

    _u = new real [(size_t) (TORD + TORD)];
    _v = _u + TORD;

  } else {
    if   (_p.z < 0.0) _p.z = Lz - fmod (fabs (_p.z), Lz);
    else              _p.z = fmod (_p.z, Lz);
    
    _u = new real [(size_t) (TORD + TORD + TORD)];
    _v = _u + TORD;
    _w = _v + TORD;
  }
}


void FluidParticle::integrate ()
// ---------------------------------------------------------------------------
// Integrate massless particle's position using predictor--corrector
// scheme.  If particles leave 2D mesh they are marked by setting E =
// 0.  For 3D, they get put back into the fundamental period of the
// solution if they leave in the z-direction.
//
// NB: The domain velocity fields must be in Fourier space prior to call (3D).
// ---------------------------------------------------------------------------
{
  if (!_E) return;
#if defined (DEBUG)
  const char routine[] = "FluidParticle::integrate";
#endif
  register integer i;
  const integer    N     = min (++_step, TORD);
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

    _u[0] = D -> u[0] -> probe (_E, _r, _s, (integer) 0);
    _v[0] = D -> u[1] -> probe (_E, _r, _s, (integer) 0);

    xp = _p.x;
    yp = _p.y;
    for (i = 0; i < N; i++) {
      xp += P_coeff[i] * _u[i];
      yp += P_coeff[i] * _v[i];
    }

    if (!_E -> locate (xp, yp, _r, _s)) {
      _E = 0;
      for (i = 0; i < NEL; i++) {
	_r = _s = 0.0;
	if (D -> elmt[i] -> locate (xp, yp, _r, _s, guess)) {
	  _E = D -> elmt[i];
	  break;
	}
      }
      if (!_E) {
#if defined (DEBUG)
	if ((int) Femlib::value("VERBOSE") > 3) {
	  char     str[StrMax];
	  sprintf (str, "Particle %1d at (%f, %f, %f) left mesh",
		   _id, _p.x, _p.y, _p.z);
	  message (routine, str, WARNING);
	}
#endif
	return;
      }
    }

    // -- Corrector.

    up = D -> u[0] -> probe (_E, _r, _s, (integer) 0);
    vp = D -> u[1] -> probe (_E, _r, _s, (integer) 0);

    _p.x += C_coeff[0] * up;
    _p.y += C_coeff[0] * vp;
    for (i = 1; i < NP; i++) {
      _p.x += C_coeff[i] * _u[i - 1];
      _p.y += C_coeff[i] * _v[i - 1];
    }

    if (!_E -> locate (_p.x, _p.y, _r, _s)) {
      _E = 0;
      for (i = 0; i < NEL; i++) {
	_r = _s = 0.0;
	if (D -> elmt[i] -> locate (_p.x, _p.y, _r, _s, guess)) {
	  _E = D -> elmt[i];
	  break;
	}
      }
      if (!_E) {
#if defined (DEBUG)
	if ((int) Femlib::value("VERBOSE") > 3) {
	  char     str[StrMax];
	  sprintf (str, "Particle %1d at (%f, %f, %f) left mesh",
		   _id, _p.x, _p.y, _p.z);
	  message (routine, str, WARNING);
	}
#endif
	return;
      }
    }

    // -- Maintain multilevel storage.

    rollv (_u, TORD);
    rollv (_v, TORD);

  } else {			// -- 3D integration.
    
    // -- Predictor.

    _u[0] = D -> u[0] -> probe (_E, _r, _s, _p.z);
    _v[0] = D -> u[1] -> probe (_E, _r, _s, _p.z);
    _w[0] = D -> u[2] -> probe (_E, _r, _s, _p.z);

    xp = _p.x;
    yp = _p.y;
    zp = _p.z;
    for (i = 0; i < N; i++) {
      xp += P_coeff[i] * _u[i];
      yp += P_coeff[i] * _v[i];
      zp += P_coeff[i] * _w[i];
    }

    if (!_E -> locate (xp, yp, _r, _s)) {
      _E = 0;
      for (i = 0; i < NEL; i++) {
	_r = _s = 0.0;
	if (D -> elmt[i] -> locate (xp, yp, _r, _s, guess)) {
	  _E = D -> elmt[i];
	  break;
	}
      }
      if (!_E) return;
    }

    // -- Corrector.

    up = D -> u[0] -> probe (_E, _r, _s, zp);
    vp = D -> u[1] -> probe (_E, _r, _s, zp);
    wp = D -> u[2] -> probe (_E, _r, _s, zp);

    _p.x += C_coeff[0] * up;
    _p.y += C_coeff[0] * vp;
    _p.z += C_coeff[0] * wp;
    for (i = 1; i < NP; i++) {
      _p.x += C_coeff[i] * _u[i - 1];
      _p.y += C_coeff[i] * _v[i - 1];
      _p.z += C_coeff[i] * _w[i - 1];
    }

    if (!_E -> locate (_p.x, _p.y, _r, _s)) {
      _E = 0;
      for (i = 0; i < NEL; i++) {
	_r = _s = 0.0;
	if (D -> elmt[i] -> locate (_p.x, _p.y, _r, _s, guess)) {
	  _E = D -> elmt[i];
	  break;
	}
      }
      if (!_E) return;
    }
    if   (_p.z < 0.0) _p.z = Lz - fmod (fabs (_p.z), Lz);
    else              _p.z = fmod (_p.z, Lz);

    // -- Maintain multilevel storage.

    rollv (_u, TORD);
    rollv (_v, TORD);
    rollv (_w, TORD);
  }
}
