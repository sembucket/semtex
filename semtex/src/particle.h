#ifndef PARTICLE_H
#define PARTICLE_H


class FluidParticle
// ===========================================================================
// Class used to locate and integrate positions of massless particles.
// ===========================================================================
{
public:
  FluidParticle (Domain*, const integer, Point&);
  void            integrate (); 
  integer         ID        () const { return _id;     }
  const  Element* inMesh    () const { return _E;      }
  const  Point&   location  () const { return _p;      } 
  static integer  IDMax     ()       { return _ID_MAX; }

private:
  integer        _id  ;		// Numeric tag.
  integer        _step;		// Number of integration steps.
  const Element* _E   ;	        // Pointer to the element particle is in.
  real           _r   ;	        // Corresponding "r" location within element.
  real           _s   ;		// likewise for "s".
  Point          _p   ;		// Physical space location.
  real*          _u   ;		// Multilevel "x" velocity storage.
  real*          _v   ;		// Multilevel "y" velocity storage.
  real*          _w   ;		// Multilevel "z" velocity storage.

  static Domain* _Dom    ;	// Velocity fields and class functions.
  static integer _NCOM   ;	// Number of velocity components.
  static integer _NEL    ;	// Number of elements in mesh.
  static integer _NZ     ;	// Number of z planes.
  static integer _TORD   ;	// Order of N--S timestepping.
  static integer _ID_MAX ;	// Highest issued id.
  static real*   _P_coeff;	// Integration (predictor) coefficients.
  static real*   _C_coeff;	// Integration (corrector) coefficients.
  static real    _DT     ;	// Time step.
  static real    _Lz     ;	// Periodic length.
};

#endif
