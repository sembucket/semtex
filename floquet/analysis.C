//////////////////////////////////////////////////////////////////////////////
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

#include <unistd.h>
#include <sem.h>


Analyser::Analyser (Domain* D   ,
		    FEML*   file) :
// ---------------------------------------------------------------------------
// Set up.
//
// Fluid particle files (session.par) are dealt with here.
// Each line is of form
// #     tag  time  x     y      z
//       1    0.0   1.0   10.0   0.5.
// Output is of the same form, called session.trk.
//
// NB: Particle tracking is broken for multiprocessor application.
//
// History points are also set up here.  They are nominated in the
// optional HISTORY section of the session file.  Output is to
// session.his.
// ---------------------------------------------------------------------------
  _src (D)
{
  const char routine[] = "Analyser::Analyser";
  char       str[StrMax];

  cout << setprecision (6);

  // -- Set up for history points: open files, create points.

  if (file -> seek ("HISTORY")) {
    integer        i, id, num = 0;
    const integer  NH = file -> attribute ("HISTORY", "NUMBER");
    const Element* E;
    HistoryPoint*  H;
    real           r, s, x, y, z;
    
    for (i = 0; i < NH; i++) {
      file -> stream() >> id >> x >> y >> z;
      if ((E = HistoryPoint::locate (x, y, D -> elmt, r, s))) {
	H = new HistoryPoint (id, E, r, s, x, y, z);
	_history.insert (_history.end(), H);
	num++;
      } else {
	sprintf (str, "History point at (%f, %f, %f) not in mesh", x, y, z);
	message (routine, str, WARNING);
      }
    }
    
    _his_strm.open (strcat (strcpy (str, _src -> name), ".his"));
    _his_strm.setf (ios::scientific, ios::floatfield);
    _his_strm.precision (6);
    if (!_his_strm) message (routine, "can't open history file", ERROR);
  }

  // -- Set up for output of modal energies every IO_CFL steps.

  _mdl_strm.open (strcat (strcpy (str, _src -> name), ".mdl"), ios::out); 
  _mdl_strm << "#     Time Mode         Energy" << endl
	   << "# ----------------------------" << endl;

  // -- Dump run information to file.

  ofstream runfile (strcat (strcpy (str, _src -> name), ".run"), ios::out);
  gethostname (str, StrMax);
  runfile << "-- Host                    : " << str << endl;
  runfile << "   PID                     : " << getpid() << endl << endl;

  D -> report (runfile);
  
  runfile.close();
}


void Analyser::analyse (AuxField** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.  If SPAWN was set, add more particles at
// original absolute positions.
// ---------------------------------------------------------------------------
{
  const integer cflstep = Femlib::ivalue ("IO_CFL");

  // -- Run information update.

  cout << "Step: " << _src -> step << "  Time: " << _src -> time << endl;

  // -- CFL, energy, divergence information.

  if (cflstep && !(_src -> step % cflstep)) {
    modalEnergy ();
    estimateCFL ();
    divergence  (work);
  }

  // -- Periodic dumps and global information.
  
  const bool periodic = !(_src -> step %  Femlib::ivalue("IO_HIS")) ||
                        !(_src -> step %  Femlib::ivalue("IO_FLD")) ;
  const bool final    =   _src -> step == Femlib::ivalue("N_STEP");
  const bool state    = periodic || final;

  if (state) {

    // -- Output history point data.
      
    register integer  i, j;
    const integer     NH = _history.size();
    const integer     NF = _src -> u.size();
    HistoryPoint*     H;
    vector<real>      tmp (NF);
    vector<AuxField*> u   (NF);

    for (i = 0; i < NF; i++)
      u[i] = _src -> u[i];

    for (i = 0; i < NH; i++) {
      H = _history[i];

      H -> extract (u, &tmp[0]);

      _his_strm << setw(4) << H->ID() << " " << setw(14) << _src->time << " ";
      for (j = 0; j < NF; j++) _his_strm << setw(15) << tmp[j];
      _his_strm << endl;
    }
  }

  _src -> dump();
}


void Analyser::modalEnergy ()
// ---------------------------------------------------------------------------
// Print out modal energies per unit area, output by root processor.
// ---------------------------------------------------------------------------
{
  const integer NC = Geometry::nPert();
  real          ek = 0.0;

  for (integer i = 0; i < NC; i++) ek += _src -> u[i] -> mode_L2 (0);

  _mdl_strm << setw(10) << _src -> time 
	    << setw( 5) << 1
	    << setw(15) << ek
	    << endl;
}


void Analyser::divergence (AuxField** Us) const
// ---------------------------------------------------------------------------
// Print out the velocity field's divergence energy per unit area.  Us
// is used as work area.
// ---------------------------------------------------------------------------
{
  const integer NC = Geometry::nPert();
  integer       i;

  if (Geometry::system() == Geometry::Cartesian) {
    for (i = 0; i < NC; i++) {
      *Us[i] = *_src -> u[i];
      Us[i] -> gradient (i);
    }
  } else {
    for (i = 0; i < NC; i++) *Us[i] = *_src -> u[i];
    Us[1] -> mulY();
    for (i = 0; i < NC; i++)  Us[i] -> gradient (i);
    Us[1] -> divY();
    if (NC == 3) Us[2] -> divY();
  }

  if (Geometry::problem() == Geometry::O2_3D_SYMM) *Us[2] *= -1.0;

  for (i = 1; i < NC; i++) *Us[0] += *Us[i];

  cout << "-- Divergence Energy: " << Us[0] -> mode_L2 (0) << endl;
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
  integer    percent;

  CFL_dt = max (_src -> u[0] -> CFL (0), _src -> u[1] -> CFL (1));
  if (Geometry::nPert() == 3) CFL_dt = max (CFL_dt, _src -> u[2] -> CFL (2));

  dt_max  = SAFETY * CFL_max / CFL_dt;
  percent = static_cast<integer>(100.0 * dt / dt_max);

  cout << "-- CFL: "     << CFL_dt * dt;
  cout << ", dt (max): " << dt_max;
  cout << ", dt (set): " << dt;
  cout << " ("           << percent << "%)" << endl;
}
