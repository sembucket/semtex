///////////////////////////////////////////////////////////////////////////////
// analysis.C: implement Analyser class for NS-type solvers.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// This deals with output of runtime information such as step numbers,
// CFL estimation, modal energies, etc.  If set, also output history
// point and particle track information.
//
// It is assumed that the first 2 or 3 (for 3D) entries in the Domain
// u vector are velocity fields.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>


Analyser::Analyser (Domain* D   ,
		    FEML*   file) :
// ---------------------------------------------------------------------------
// Set up.
// ---------------------------------------------------------------------------
  _src (D)
{
  char str[StrMax];

  cout << setprecision (6);

  // -- Initialize averaging.

  if (Femlib::ivalue ("AVERAGE")) {
    _stats = new Statistics (D);
  } else                              
    _stats = 0;

  // -- Set up for output of modal energies every IO_CFL steps if 3D.

  if (Geometry::nDim() == 3) {
    strcat (strcpy (str, _src -> name), ".mdl");
    ROOTONLY {
      _mdl_strm.open (str, ios::out); 
      _mdl_strm <<
"#         Time          Mode0          ModeC       ModeC.Re       ModeC.Im"
<< endl
	       <<
"# ------------------------------------------------------------------------" 
<< endl;
      _mdl_strm.setf (ios::scientific, ios::floatfield);
      _mdl_strm.precision (8);
    }
  }
}


void Analyser::analyse (AuxField** work0,
			AuxField** work1)
// ---------------------------------------------------------------------------
// Step-by-step processing.  If SPAWN was set, add more particles at
// original absolute positions.
// ---------------------------------------------------------------------------
{
  const int_t verbose = Femlib::ivalue ("VERBOSE");
  const int_t cflstep = Femlib::ivalue ("IO_CFL");

  // -- Step-by-step updates.

  cout << "Step: " << _src -> step << "  Time: " << _src -> time << endl;

  // -- CFL, energy, divergence information.

  if (cflstep && !(_src -> step % cflstep)) {
    this -> modalEnergy ();
    this -> estimateCFL ();
    this -> divergence  (work0);
  }

  // -- Periodic dumps and global information.
  
  const bool periodic = !(_src->step %  Femlib::ivalue("IO_HIS")) ||
                        !(_src->step %  Femlib::ivalue("IO_FLD"));
  const bool final    =   _src->step == Femlib::ivalue("N_STEP");
  const bool state    = periodic || final;

  if (state) {
     
    // -- Statistical analysis.

    if (_stats) _stats -> update (work0, work1);
  }

  // -- Field and statistical dumps.

  _src -> dump ();
  if (_stats) {
    char filename[StrMax];
    strcat (strcpy (filename, _src -> name), ".avg");
    _stats -> dump (filename);
  }
}


void Analyser::modalEnergy ()
// ---------------------------------------------------------------------------
// Print out modal energies per unit area, output by root processor.
// ---------------------------------------------------------------------------
{
  const int_t    DIM   = Geometry::nDim();
  register int_t i;
  real_t         re, im, ek[4];

  Veclib::zero (4, ek, 1);

  for (i = 0; i < DIM; i++) {
    _src -> u[i] -> mode_en (0, re, im);
    ek[0] += re;
    _src -> u[i] -> mode_en (1, re, im);
    ek[1] += re + im;
    ek[2] += re;
    ek[3] += im;
  }

  _mdl_strm << setw(10) << _src -> time 
	    << setw(15) << ek[0]
	    << setw(15) << ek[1]
	    << setw(15) << ek[2]
	    << setw(15) << ek[3]
	    << endl;
}


void Analyser::divergence (AuxField** Us) const
// ---------------------------------------------------------------------------
// Print out the velocity field's divergence energy per unit area.  Us
// is used as work area.
// ---------------------------------------------------------------------------
{
  const Geometry::CoordSys space = Geometry::system();

  const int_t    DIM = Geometry::nDim();
  const int_t    N   = Geometry::nModeProc();
  const real_t   Lz  = Femlib::value ("TWOPI / BETA");
  register int_t i, m;
  real_t         L2 = 0.0;

  if (space == Geometry::Cartesian) {

    for (i = 0; i < DIM; i++) {
      *Us[i] = *_src -> u[i];
      Us[i] -> gradient (i);
    }

  } else {			// -- Cylindrical.

    for (i = 0; i < DIM; i++) *Us[i] = *_src -> u[i];
    Us[1] -> mulY();
    for (i = 0; i < DIM; i++)  Us[i] -> gradient (i);
    Us[1] -> divY();
    if (DIM == 3) Us[2] -> divY();

  }

  for (i = 1; i < DIM; i++) *Us[0] += *Us[i];

  for (m = 0; m < N; m++) L2 += Us[0] -> mode_L2 (m);

  L2 /= Lz;

  cout << "-- Divergence Energy: " << L2 << endl;
}


void Analyser::estimateCFL () const
// ---------------------------------------------------------------------------
// Estimate and print the peak CFL number, based on zero-mode velocities.
// ---------------------------------------------------------------------------
{
  const real_t CFL_max = 0.7;	// -- Approximate maximum for scheme.
  const real_t SAFETY  = 0.9;	// -- Safety factor.
  const real_t dt      = Femlib::value ("D_T");
  real_t       CFL_dt, dt_max;
  int_t        percent;

  CFL_dt = max (_src -> u[0] -> CFL (0), _src -> u[1] -> CFL (1));
  if (Geometry::nDim() == 3) CFL_dt = max (CFL_dt, _src -> u[2] -> CFL (2));

  dt_max  = SAFETY * CFL_max / CFL_dt;
  percent = static_cast<int_t>(100.0 * dt / dt_max);

  cout << "-- CFL: "     << CFL_dt * dt;
  cout << ", dt (max): " << dt_max;
  cout << ", dt (set): " << dt;
  cout << " ("           << percent << "%)" << endl;
}
