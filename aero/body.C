//////////////////////////////////////////////////////////////////////////////
// body.C: implement Body class for aeroelastic coupling.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <aero.h>

#if defined(__DECCXX)
  #pragma define_template roll<real>
  #pragma define_template min<integer>
  #pragma define_template max<integer>
#endif

static AxisMotion* createAxis (char*);


Body::Body (const char* session)
// ---------------------------------------------------------------------------
// Input file sets up a body.  If no file, there is no motion on either
// axis, but forces exerted on walls are computed and printed.
//
// Input file format: {entry not allowed in file}/(treated as a comment).
// -- BOF --
// Title string.                          {Just a header, not used.}
// Another string or blank line.          (Cosmetic padding.)
// x-axis cosine     ampl freq phase      (3 real values in ASCII.)
// y-axis feedback   mass freq zeta       (3 real values in ASCII.)
// x-state           pos  vel             (2 real values, default to zero.)
// x-state           pos  vel             (Repetitions allowed, last used.)
// y-state           pos  vel             (X & Y states can be intermixed.)
// -- EOF --
//
// Notes:
// 1. The x-axis and y-axis lines must both appear, in order.  The
//    two first (title) lines must also appear.  Lower case only.
// 2. To hold an axis stationary use "fixed".
// 3. Motion can be prescribed using type "function", in which case 3
//    function strings defining position, velocity and acceleration follow.
// 4. State information is ignored except in case of "feedback" motion.
//    Only the last line of state information is used for each axis.
// 5. Timebase for prescribed and sinusoidal motions comes from restart file
//    and prescribed time step.
// 6. For cosinusoidal motions:
//    Amplitude is zero-to-peak amplitude.
//    Frequency freq is given in dimensionless frequency units.
//    Phase angle is given in radians.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Body::Body";
  char       s[StrMax], err[StrMax];
  ifstream   file (strcat (strcpy (s, session), ".bdy"));

  if (file.fail ()) {		// -- Default action for no ".bdy" file.
    axis[0] = createAxis ("fixed");
    axis[1] = createAxis ("fixed");

  } else {			// -- File exists.  Create 2 axes.
    file.getline (s, StrMax);
    file.getline (s, StrMax);	// -- Strip header.

    char* tok;
    char  sep[] = " \t";

    file.getline (s, StrMax);
    tok = strtok (s, sep);
    if (strstr (tok, "x-axis")) {
      tok = strtok (0, "\0");
      axis[0] = createAxis (tok);
    } else {
      ostrstream (err, StrMax) 
	<< "expected x-axis specifier, got: " << s << ends;
      message (routine, err, ERROR);
    }

    file.getline (s, StrMax);
    tok = strtok (s, sep);
    if (strstr (tok, "y-axis")) {
      tok = strtok (0, "\0");
      axis[1] = createAxis (tok);
    } else {
      ostrstream (err, StrMax) 
	<< "expected y-axis specifier, got: " << s << ends;
      message (routine, err, ERROR);
    }
  }
  
  state.x = state.y = state.z = 0.0;
  traction[0] = traction[1] = state;

  // -- Take special action to install motion state variables for SMD axis.

  axis[0] -> describe (s);
  if (strstr (s, "mass")) ((SMD*) axis[0]) -> setState (file, "x-state");

  axis[1] -> describe (s);
  if (strstr (s, "mass")) ((SMD*) axis[1]) -> setState (file, "y-state");

  // -- Echo input to cout.

  cout << "-- Body axis information:" << endl;
  axis[0] -> describe (s);
  cout << "   x axis: " << s << endl;
  axis[1] -> describe (s);
  cout << "   y axis: " << s << endl;

  file.close();
}


Vector Body::acceleration ()
// ---------------------------------------------------------------------------
// Return acceleration on both axes.
// ---------------------------------------------------------------------------
{
  state.x = axis[0] -> getA ();
  state.y = axis[1] -> getA ();

  return state;
}


Vector Body::velocity ()
// ---------------------------------------------------------------------------
// Return velocity on both axes.
// ---------------------------------------------------------------------------
{
  state.x = axis[0] -> getV ();
  state.y = axis[1] -> getV ();

  return state;
}


Vector Body::position ()
// ---------------------------------------------------------------------------
// Return position on both axes.
// ---------------------------------------------------------------------------
{
  state.x = axis[0] -> getX ();
  state.y = axis[1] -> getX ();

  return state;
}


void Body::move (const integer step)
// ---------------------------------------------------------------------------
// Update estimates of position, velocity and acceleration at next time level.
// Note that for SMD class, acceleration is updated by Body::force instead.
// ---------------------------------------------------------------------------
{
  axis[0] -> move (step);
  axis[1] -> move (step);
}


Vector Body::force (const Domain* D)
// ---------------------------------------------------------------------------
// Compute and install pressure, viscous and total forces exerted on body.
// Return total force.
// ---------------------------------------------------------------------------
{
  const integer DIM = Geometry::nDim();

  traction[0] = Field::normalTraction  (D->u[DIM]);
  traction[1] = Field::tangentTraction (D->u[0], D->u[1]);

  traction[2].x = traction[0].x + traction[1].x;
  traction[2].y = traction[0].y + traction[1].y;

  axis[0] -> force (traction[2].x);
  axis[1] -> force (traction[2].y);

  return traction[2];
}


ostream& operator << (ostream& S, Body* B)
// ---------------------------------------------------------------------------
// Print up motion state variables (x, xdot, xddot, Fvis, Fpre, Ftot).
// ---------------------------------------------------------------------------
{
  char   s[StrMax];
  Vector x, xdot, xddot;

  x     = B->position     ();
  xdot  = B->velocity     ();
  xddot = B->acceleration ();

  sprintf (s, " %#10.6g %#10.6g %#10.6g %#10.6g %#10.6g %#10.6g"
              " %#10.6g %#10.6g %#10.6g %#10.6g %#10.6g %#10.6g",
	               x.x,          xdot.x,         xddot.x,
	   B->traction[0].x, B->traction[1].x, B->traction[2].x,
	               x.y,          xdot.y,         xddot.y, 
	   B->traction[0].y, B->traction[1].y, B->traction[2].y);
  S << s;

  return S;
}


Cosine::Cosine (char* s)
// ---------------------------------------------------------------------------
// Copy strings that, when interpreted, return amplitude, frequency and phase
// for cosinusiodal motion.
// 
// Amplitude is zero-to-peak amplitude.
// Frequency freq is given in dimensionless frequency units.
// Phase angle is given in radians.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Cosine::Cosine";
  char       err[StrMax];
  char*      tok;
  char       sep[] = " \t";

  tok = strtok (s, sep); strcpy (amplitude, tok);
  tok = strtok (0, sep); strcpy (frequency, tok);
  tok = strtok (0, sep); strcpy (phaseangl, tok);

  if (   (strlen (amplitude) == 0)
      || (strlen (frequency) == 0)
      || (strlen (phaseangl) == 0)) {
    ostrstream (err, StrMax)
      << "couldn't parse 3 real values from string: " << s << ends;
    message (routine, err, ERROR);
  }
  
  move ();
}


void Cosine::move (const integer dummi)
// ---------------------------------------------------------------------------
// Update state variables.  Dummy input variable not used.
// ---------------------------------------------------------------------------
{
  real mag   = Femlib::value (amplitude);
  real omega = Femlib::value ("TWOPI") * Femlib::value (frequency);
  real angle = omega * Femlib::value ("t") + Femlib::value (phaseangl);

  pos =                  mag * cos (angle);
  vel =         -omega * mag * sin (angle);
  acc = -omega * omega * mag * cos (angle);
}

 
void Cosine::describe (char* s)
// ---------------------------------------------------------------------------
// Return description in s.
// ---------------------------------------------------------------------------
{
  ostrstream (s, StrMax)
    << "cosine: "
      << amplitude 
	<<" * cos (2PI * " 
	  << frequency
	    << " * t + " 
	      << phaseangl
		<< ")"
		  << ends;
}


Function::Function (char* s)
// ---------------------------------------------------------------------------
// String s contains three function strings; in order: pos, vel, acc.
// Tokenize them and do first parse to set internal state variables.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Function::Function";
  char       err[StrMax];
  char*      tok;
  char       sep[] = " \t";

  tok = strtok (s, sep); strcpy (position,     tok);
  tok = strtok (0, sep); strcpy (velocity,     tok);
  tok = strtok (0, sep); strcpy (acceleration, tok);

  if (   (strlen (position)     == 0)
      || (strlen (velocity)     == 0)
      || (strlen (acceleration) == 0)) {
    ostrstream (err, StrMax)
      << "couldn't parse 3 real values from string: " << s << ends;
    message (routine, err, ERROR);
  }

  pos = Femlib::value (position);
  vel = Femlib::value (velocity);
  acc = Femlib::value (acceleration);
}


void Function::move (const integer dummi)
// ---------------------------------------------------------------------------
// Update state variables.  Neither input is used.
// ---------------------------------------------------------------------------
{
  pos = Femlib::value (position);
  vel = Femlib::value (velocity);
  acc = Femlib::value (acceleration);
}

 
void Function::describe (char* s)
// ---------------------------------------------------------------------------
// Return description in s.
// ---------------------------------------------------------------------------
{
  ostrstream (s, StrMax)
    << "function:"
      << " pos: "
	<< position
	  << " vel: "
	    << velocity
	      << " acc: "
		<< acceleration
		  << ends;
}


SMD::SMD (char* s)
// ---------------------------------------------------------------------------
// String s contains three function strings; in order: mass per unit length,
// natural frequency, and damping ratio.
// Tokenize them.  Initialize state motion variables.
// ---------------------------------------------------------------------------
{
  const char routine[] = "SMD::SMD";
  char       err[StrMax], sep[] = " \t";;
  char*      tok;

  tok = strtok (s, sep); strcpy (mass, tok);
  tok = strtok (0, sep); strcpy (natf, tok);
  tok = strtok (0, sep); strcpy (zeta, tok);

  if ((strlen (mass) == 0) || (strlen (natf) == 0) || (strlen (zeta) == 0)) {
    ostrstream (err, StrMax)
      << "couldn't parse 3 real values from string: " << s << ends;
    message (routine, err, ERROR);
  }

  x    = new real [Integration::OrderMax];
  xdot = new real [Integration::OrderMax];
  f    = new real [Integration::OrderMax];

  pos = vel = acc = x[0] = xdot[0] = f[0] = 0.0;
}


void SMD::move (const integer step)
// ---------------------------------------------------------------------------
// Update state variables using stiffly-stable time integration scheme.
// 
// On exit, position and velocity have been integrated to end of new step.
// ---------------------------------------------------------------------------
{
  const real    w    = Femlib::value ("TWOPI") * Femlib::value (natf);
  const real    z    = Femlib::value (zeta);
  const real    dt   = Femlib::value ("D_T");
  const integer Jmax = (integer) Femlib::value ("N_TIME");

  integer       q, Je = min (max ((integer) 1, step), Jmax);
  real          *alpha, *betaDt;

  vector<real> work (2 * Je + 1);

  alpha  = work();
  betaDt = alpha + Je + 1;

  Integration::StifflyStable (Je, alpha);
  Integration::Extrapolation (Je, betaDt);
  Veclib::smul (Je, dt, betaDt, 1, betaDt, 1);

  // -- Integrate motion state variables.

  vel = 0.0;
  for (q = 0; q < Je; q++) {
    vel -= alpha[q + 1] *     xdot[q];
    vel += betaDt[q] *           f[q];
    vel -= betaDt[q] * sqr (w) * x[q];
  }
  vel /= alpha[0] + 2.0 * z * w * dt;

  pos = dt * vel;
  for (q = 0; q < Je; q++)
    pos -= alpha[q+1] * x[q];
  pos /= alpha[0];

  // -- Maintain motion state variable FIFO storage.
  
  rollv (x,    Jmax);    x[0] = pos;
  rollv (xdot, Jmax); xdot[0] = vel;
}


void SMD::force (const real F)
// ---------------------------------------------------------------------------
// Update force per unit mass extrapolation vector, estimate acceleration
// at current time level.
// 
// F is the tractive force on this axis at end of current time step.
// ---------------------------------------------------------------------------
{
  const real    m    = Femlib::value (mass);
  const real    w    = Femlib::value ("TWOPI") * Femlib::value (natf);
  const real    z    = Femlib::value (zeta);
  const integer Jmax = (integer) Femlib::value ("N_TIME");

  rollv (f, Jmax);
  f[0] = F / m;

  acc = f[0] - 2.0 * z * w * vel - sqr (w) * pos;
}

 
void SMD::describe (char* s)
// ---------------------------------------------------------------------------
// Return description in s.
// ---------------------------------------------------------------------------
{
  ostrstream (s, StrMax)
    << "feedback:"
      << " mass: "
	<< mass
	  << " freq: "
	    << natf
	      << " zeta: "
		<< zeta
		  << " pos: "
		    << x[0]
		      << " vel: "
			<< xdot[0]
			  << ends;
}


void SMD::setState (ifstream& file, const char* tag)
// ---------------------------------------------------------------------------
// Traverse file looking for lines beginning with tag.  Expect x, and xdot
// state variable values to follow.  Install last corresponding line in file.
// ---------------------------------------------------------------------------
{
  const char routine[] = "SMD::setState";
  char       s[StrMax], err[StrMax], sep[] = " \t";
  char*      head;
  char*      tail;

  file.clear ();
  file.seekg (0);

  while (file.getline (s, StrMax)) {
    head = strtok (s, sep);
    tail = strtok (0, "\0");
    if (head && strstr (head, tag)) {
      istrstream t (tail, strlen (tail));
      t >> x[0] >> xdot[0];
      if (t.bad ()) {
	ostrstream (err, StrMax)
	  << "couldn't parse position & velocity for "
	    << tag 
	      << " from: " 
		<< s 
		  << ends;
	message (routine, err, ERROR);
      }
    }
  }
}


AxisMotion* createAxis (char* s)
// ---------------------------------------------------------------------------
// Factory for different axis-motion types.  Return pointer to abstract
// base class.  
// ---------------------------------------------------------------------------
{
  AxisMotion* base = 0;

  const char  routine[]  = "createAxis";
  char        err[StrMax], buf[StrMax], sep[] = " \t";
  char*       kind;
  char*       tail;
  
  strcpy (buf, s);
  kind = strtok (buf, sep);
  tail = strtok (0, "\0");

  if (strstr (kind, "fixed"))
    base = new Fixed;

  else if (strstr (kind, "cosine"))
    base = new Cosine (tail);

  else if (strstr (kind, "function"))
    base = new Function (tail);

  else if (strstr (kind, "feedback"))
    base = new SMD (tail);

  else {
    ostrstream (err, StrMax)
      << "couldn't parse a known axis kind from string: " << s << ends;
    message (routine, err, ERROR);
  }

  return base;
}
