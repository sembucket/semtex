///////////////////////////////////////////////////////////////////////////////
// analysis.C: implement Analyser class for NS-type solvers.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
//
// This deals with output of runtime information such as step numbers,
// CFL estimation, modal energies, etc.  If set, also output history
// point and particle track information.
//
// It is assumed that the first 2 or 3 (for 3D) entries in the Domain
// u vector are velocity fields.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


Analyser::Analyser (Domain* D   ,
		    FEML*   file) :
// ---------------------------------------------------------------------------
// Set up.
// ---------------------------------------------------------------------------
  src (D)
{
  char str[StrMax];

  cout << setprecision (6);

  // -- Initialize averaging.

  if ((integer) Femlib::value ("AVERAGE")) {
    vector<AuxField*> extra (0);
    stats = new Statistics (D, extra);
  } else                              
    stats = 0;

  // -- Set up for output of modal energies every IO_CFL steps if 3D.

  if (Geometry::nDim() == 3) {
    strcat (strcpy (str, src -> name), ".mdl");
    ROOTONLY {
      mdl_strm.open (str, ios::out); 
      mdl_strm <<
"#         Time          Mode0          ModeC       ModeC.Re       ModeC.Im"
<< endl
	       <<
"# ------------------------------------------------------------------------" 
<< endl;
      mdl_strm.setf (ios::scientific, ios::floatfield);
      mdl_strm.precision (8);
    }
  }
}


void Analyser::analyse (AuxField** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.  If SPAWN was set, add more particles at
// original absolute positions.
// ---------------------------------------------------------------------------
{
  const integer verbose = (integer) Femlib::value ("VERBOSE");
  const integer cflstep = (integer) Femlib::value ("IO_CFL");

  // -- Step-by-step updates.

  cout << "Step: " << src -> step << "  Time: " << src -> time << endl;

  // -- CFL, energy, divergence information.

  if (cflstep && !(src -> step % cflstep)) {
    modalEnergy ();
    estimateCFL ();
    divergence  (work);
  }

  // -- Periodic dumps and global information.
  
  const integer periodic = !(src->step %  (integer) Femlib::value("IO_HIS")) ||
                           !(src->step %  (integer) Femlib::value("IO_FLD"));
  const integer final    =   src->step == (integer) Femlib::value("N_STEP");
  const integer state    = periodic || final;

  if (state) {
     
    // -- Statistical analysis.

    if (stats) stats -> update (work);
  }

  // -- Field and statistical dumps.

  src -> dump ();
  if (stats) stats -> dump();
}


void Analyser::modalEnergy ()
// ---------------------------------------------------------------------------
// Print out modal energies per unit area, output by root processor.
// ---------------------------------------------------------------------------
{
  const integer    DIM   = Geometry::nDim();
  register integer i;
  real             re, im, ek[4];

  Veclib::zero (4, ek, 1);

  for (i = 0; i < DIM; i++) {
    src -> u[i] -> mode_en (0, re, im);
    ek[0] += re;
    src -> u[i] -> mode_en (1, re, im);
    ek[1] += re + im;
    ek[2] += re;
    ek[3] += im;
  }

  mdl_strm << setw(10) << src -> time 
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

  const integer    DIM = Geometry::nDim();
  const integer    N   = Geometry::nModeProc();
  const real       Lz  = Femlib::value ("TWOPI / BETA");
  register integer i, m;
  real             L2 = 0.0;

  if (space == Geometry::Cartesian) {

    for (i = 0; i < DIM; i++) {
      *Us[i] = *src -> u[i];
      Us[i] -> gradient (i);
    }

  } else {			// -- Cylindrical.

    for (i = 0; i < DIM; i++) *Us[i] = *src -> u[i];
    Us[1] -> mulR();
    for (i = 0; i < DIM; i++)  Us[i] -> gradient (i);
    Us[1] -> divR();
    if (DIM == 3) Us[2] -> divR();

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
  const real CFL_max = 0.7;	// -- Approximate maximum for scheme.
  const real SAFETY  = 0.9;	// -- Saftey factor.
  const real dt      = Femlib::value ("D_T");
  real       CFL_dt, dt_max;
  int        percent;

  CFL_dt = max (src -> u[0] -> CFL (0), src -> u[1] -> CFL (1));
  if (Geometry::nDim() == 3) CFL_dt = max (CFL_dt, src -> u[2] -> CFL (2));

  dt_max  = SAFETY * CFL_max / CFL_dt;
  percent = (int) (100.0 * dt / dt_max);

  cout << "-- CFL: "     << CFL_dt * dt;
  cout << ", dt (max): " << dt_max;
  cout << ", dt (set): " << dt;
  cout << " ("           << percent << "%)" << endl;
}
