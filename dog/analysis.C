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
#include <unistd.h>

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
  src (D)
{
  const char routine[] = "Analyser::Analyser";
  char       str[StrMax];

  cout << setprecision (6);

  // -- Set up for history points: open files, create points.

  if (file -> seek ("HISTORY")) {
    integer              i, id, num = 0;
    const integer        NH = file -> attribute ("HISTORY", "NUMBER");
    const Element*       E;
    HistoryPoint*        H;
    Stack<HistoryPoint*> stack;
    real                 r, s, x, y, z;
    
    for (i = 0; i < NH; i++) {
      file -> stream() >> id >> x >> y >> z;
      if ((E = HistoryPoint::locate (x, y, D -> elmt, r, s))) {
	H = new HistoryPoint (id, E, r, s, z);
	stack.push (H);
	num++;
      } else {
	sprintf (str, "History point at (%f, %f, %f) not in mesh", x, y, z);
	message (routine, str, WARNING);
      }
    }
    
    history.setSize (num);
    while (num--) history[num] = stack.pop();

    his_strm.open (strcat (strcpy (str, src -> name), ".his"));
    his_strm.setf (ios::scientific, ios::floatfield);
    his_strm.precision (6);
    if (!his_strm) message (routine, "can't open history file", ERROR);
  }

  // -- Initialize averaging.

  if ((integer) Femlib::value ("AVERAGE")) {
    vector<AuxField*> extra (0);
    stats = new Statistics (D, extra);
  } else                              
    stats = 0;

  // -- Set up for output of modal energies every IO_CFL steps.

  mdl_strm.open (strcat (strcpy (str, src -> name), ".mdl"), ios::out); 
  mdl_strm << "#     Time Mode         Energy" << endl
	   << "# ----------------------------" << endl;

  // -- Dump run information to file.

  ofstream runfile (strcat (strcpy (str, src -> name), ".run"), ios::out);
  int pid = getpid(); gethostname (str, StrMax);
  runfile << "host = " << str << endl << "pid  = " << pid << endl;
  runfile.close();
}


void Analyser::analyse (AuxField** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.  If SPAWN was set, add more particles at
// original absolute positions.
// ---------------------------------------------------------------------------
{
  const integer cflstep = (integer) Femlib::value ("IO_CFL");

  // -- Run information update.

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

    // -- Output history point data.
      
    register integer  i, j;
    const integer     NH = history.getSize();
    const integer     NF = src -> u.getSize();
    HistoryPoint*     H;
    vector<real>      tmp (NF);
    vector<AuxField*> u   (NF);

    for (i = 0; i < NF; i++)
      u[i] = src -> u[i];

    for (i = 0; i < NH; i++) {
      H = history[i];

      H -> extract (u, tmp());

      his_strm << setw(4) << H->ID() << " " << setw(14) << src->time << " ";
      for (j = 0; j < NF; j++) his_strm << setw(15) << tmp[j];
      his_strm << endl;
    }

    // -- Statistical analysis.

    if (stats) stats -> update (work);
  }

  // -- Field and statistical dumps.

  src -> dump();
  if (stats) stats -> dump();
}


void Analyser::modalEnergy ()
// ---------------------------------------------------------------------------
// Print out modal energies per unit area, output by root processor.
// ---------------------------------------------------------------------------
{
  const integer    DIM   = Geometry::nDim();
  const integer    N     = Geometry::nModeProc();
  const integer    nProc = Geometry::nProc();
  register integer i, m;
  vector<real>     ek (N);

  for (m = 0; m < N; m++) {
    ek[m] = 0.0;
    for (i = 0; i < DIM; i++) ek[m] += src -> u[i] -> mode_L2 (m);
  }

  for (m = 0; m < N; m++)
    mdl_strm << setw(10) << src -> time 
	     << setw( 5) << m 
	     << setw(15) << ek[m]
	     << endl;
}


void Analyser::divergence (AuxField** Us) const
// ---------------------------------------------------------------------------
// Print out the velocity field's divergence energy per unit area.  Us
// is used as work area.
// ---------------------------------------------------------------------------
{
  const Geometry::CoordSys space = Geometry::system();

  const integer    DIM = Geometry::nPert();
  const integer    N   = Geometry::nModeProc();
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

  if (DIM == 3) *Us[2] *= -1.0;	// -- Third velocity component is imaginary.

  for (i = 1; i < DIM; i++) *Us[0] += *Us[i];

  for (m = 0; m < N; m++) L2 += Us[0] -> mode_L2 (m);

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
