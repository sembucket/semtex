#ifndef AERO_H
#define AERO_H
//////////////////////////////////////////////////////////////////////////////
// aero.h: header file for non-inertial NS solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "sem.h"

class Body;


class AeroAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  AeroAnalyser (Domain*, FEML*, Body*);
  void analyse (AuxField**);

private:
  Body*    body;
  ofstream sta_strm;
  ofstream for_strm;
  void     forceDist();
};


class AxisMotion
// ===========================================================================
// Interface base class for body movement.
// ===========================================================================
{
public:
  virtual ~AxisMotion() { };

  virtual real_t getA () = 0;
  virtual real_t getV () = 0;
  virtual real_t getX () = 0;
  
  virtual void move  (const int_t)  = 0;
  virtual void force (const real_t) = 0;

  virtual void describe (char*) = 0;
};


class Fixed : public AxisMotion
// ===========================================================================
// This axis is fixed.
// ===========================================================================
{
public:
  Fixed () : zero (0) { }

  virtual real_t getA () { return zero; }
  virtual real_t getV () { return zero; }
  virtual real_t getX () { return zero; }
  
  virtual void move  (const int_t = 0) { }
  virtual void force (const real_t)    { }

  virtual void describe (char* s) { strcpy (s, "fixed"); }

 private:
  real_t zero;
};


class Cosine : public AxisMotion
// ===========================================================================
// Cosinusoidal motion x = amplitude * cos (2Pi * frequency * t + phaseangl).
// Private motion state variables are updated by calling move (and are
// initially set by constructor).
// ===========================================================================
{
public:
  Cosine (char*);

  virtual real_t getA () { return acc; }
  virtual real_t getV () { return vel; }
  virtual real_t getX () { return pos; }

  virtual void move  (const int_t = 0);
  virtual void force (const real_t) { }

  virtual void describe (char*);

private:
  char amplitude[StrMax];
  char frequency[StrMax];
  char phaseangl[StrMax];

  real_t pos;
  real_t vel;
  real_t acc;
};


class Function : public AxisMotion
// ===========================================================================
// Prescribed motion.  Constructor string has position, velocity and
// acceleration functions, for interpretation.
// Private motion state variables are updated by calling move (and are
// initially set by constructor).
// ===========================================================================
{
public:
  Function (char*);

  virtual real_t getA () { return acc; }
  virtual real_t getV () { return vel; }
  virtual real_t getX () { return pos; }

  virtual void move  (const int_t = 0);
  virtual void force (const real_t) { }

  virtual void describe (char*);

private:
  char acceleration[StrMax];
  char velocity[StrMax];
  char position[StrMax];
  
  real_t pos;
  real_t vel;
  real_t acc;
};


class SMD : public AxisMotion
// ===========================================================================
// Spring-Mass-Damper class for 1D motion.  Constructor string supplies
// mass per unit length, natural frequency and damping ratio as strings
// for interpretation.  Time integration is by "stiffly-stable" scheme.
// ===========================================================================
{
public:
  SMD (char*);

  virtual real_t getA () { return acc; }
  virtual real_t getV () { return vel; }
  virtual real_t getX () { return pos; }

  virtual void move  (const int_t );
  virtual void force (const real_t);

  virtual void describe (char*);

  void    setState (ifstream&, const char*);

private:
  char mass[StrMax];
  char natf[StrMax];
  char zeta[StrMax];
  
  real_t pos;
  real_t vel;
  real_t acc;

  real_t* xdot;			// -- Velocity state variables for integration.
  real_t* x;			// -- Position state variables.
  real_t* f;			// -- Force per unit mass.
};
  

class Body
// ===========================================================================
// Implement body class, with functions to enable coupling to Navier--Stokes
// solver.
// ===========================================================================
{
friend ostream& operator << (ostream&, Body*);
public:
  Body  (const char*);
  ~Body () { }

  Vector  acceleration ();
  Vector  velocity     ();
  Vector  position     ();
  Vector  force        (const Domain*);

  void    move         (const int_t);

private:
  AxisMotion*  axis[2];		// -- Array of motion interface classes.
  Vector       state;		// -- Dummy for output.
  Vector       traction[3];	// -- Pressure, viscous, total forces.
};

void NavierStokes  (Domain*, Body*, AeroAnalyser*);
void eddyViscosity (const Domain*, AuxField**, AuxField**, AuxField*);
#endif
