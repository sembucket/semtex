//////////////////////////////////////////////////////////////////////////////
// aero.h: header file for non-inertial NS solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


class Body;


class Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  Analyser  (Domain&, Body&);
  ~Analyser () { }

  void analyse ();

private:
  Domain&   src;
  ofstream  fld_strm;

  Body&     body;
  ofstream  sta_strm;

  void modalEnergy () const;
};


class AxisMotion
// ===========================================================================
// Interface base class for body movement.
// ===========================================================================
{
public:
  virtual ~AxisMotion() { };

  virtual real getA () = 0;
  virtual real getV () = 0;
  virtual real getX () = 0;
  
  virtual void move  (const int)  = 0;
  virtual void force (const real) = 0;

  virtual void describe (char*) = 0;
};


class Fixed : public AxisMotion
// ===========================================================================
// This axis is fixed.
// ===========================================================================
{
public:
  Fixed () : zero (0) { }

  virtual real getA () { return zero; }
  virtual real getV () { return zero; }
  virtual real getX () { return zero; }
  
  virtual void move  (const int = 0) { }
  virtual void force (const real)    { }

  virtual void describe (char* s) { strcpy (s, "fixed"); }

 private:
  real zero;
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

  virtual real getA () { return acc; }
  virtual real getV () { return vel; }
  virtual real getX () { return pos; }

  virtual void move  (const int = 0);
  virtual void force (const real) { }

  virtual void describe (char*);

private:
  char amplitude[StrMax];
  char frequency[StrMax];
  char phaseangl[StrMax];

  real pos;
  real vel;
  real acc;
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

  virtual real getA () { return acc; }
  virtual real getV () { return vel; }
  virtual real getX () { return pos; }

  virtual void move  (const int = 0);
  virtual void force (const real) { }

  virtual void describe (char*);

private:
  char acceleration[StrMax];
  char velocity[StrMax];
  char position[StrMax];
  
  real pos;
  real vel;
  real acc;
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

  virtual real getA () { return acc; }
  virtual real getV () { return vel; }
  virtual real getX () { return pos; }

  virtual void move  (const int );
  virtual void force (const real);

  virtual void describe (char*);

  void    setState (ifstream&, const char*);

private:
  char mass[StrMax];
  char natf[StrMax];
  char zeta[StrMax];
  
  real pos;
  real vel;
  real acc;

  real* xdot;			// -- Velocity state variables for integration.
  real* x;			// -- Position state variables.
  real* f;			// -- Force per unit mass.
};
  

class Body
// ===========================================================================
// Implement body class, with functions to enable coupling to Navier--Stokes
// solver.
// ===========================================================================
{
friend ostream& operator << (ostream&, Body&);
public:
  Body  (const char*);
  ~Body () { }

  Vector  acceleration ();
  Vector  velocity     ();
  Vector  position     ();
  Vector  force        (const Domain&);

  void    move  (const int);

private:
  AxisMotion*  axis[2];		// -- Array of motion interface classes.
  Vector       state;		// -- Dummy for output.
  Vector       traction[3];	// -- Pressure, viscous, total forces.
};

