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

  // -- Set up for particle tracking.

  ROOTONLY {
    
    // -- Open particle track file, create particles.

    ifstream pfile (strcat (strcpy (str, src -> name), ".par"));  

   if (!pfile.fail()) {
      const integer  add = (integer) Femlib::value ("SPAWN");
      integer        id;
      Point          P, *I;
      FluidParticle* F;

      par_strm.open (strcat (strcpy (str, src -> name), ".trk"));
      par_strm.setf (ios::scientific, ios::floatfield);
      par_strm.precision (6);

      while (pfile >> id >> P.x >> P.x >> P.y >> P.z) {
	F = new FluidParticle (src, id, P);
	if (!(F -> inMesh())) {
	  sprintf (str, "Particle at (%f, %f, %f) not in mesh", P.x, P.y, P.z);
	  message (routine, str, WARNING);
	} else
	  particle.add (F);
	if (add) {
	  I = new Point;
	  I -> x = P.x; I -> y = P.y; I -> z = P.z;
	  initial.add (I);
	}
      }
    }
  }

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
      if (E = HistoryPoint::locate (x, y, D -> elmt, r, s)) {
	H = new HistoryPoint (id, E, r, s, x, y, z);
	stack.push (H);
	num++;
      } else {
	sprintf (str, "History point at (%f, %f, %f) not in mesh", x, y, z);
	message (routine, str, WARNING);
      }
    }
    
    history.setSize (num);
    while (num--) history[num] = stack.pop();

    ROOTONLY {
      his_strm.open (strcat (strcpy (str, src -> name), ".his"));
      his_strm.setf (ios::scientific, ios::floatfield);
      his_strm.precision (6);
      if (!his_strm) message (routine, "can't open history file", ERROR);
    }
  }

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
      mdl_strm << "#     Time Mode         Energy" << endl
	       << "# ----------------------------" << endl;
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
  const integer add     = (integer) Femlib::value ("SPAWN") &&
                         !(src -> step % (integer) Femlib::value ("SPAWN"));

  ListIterator<FluidParticle*> p (particle);

  // -- Step-by-step updates.

  ROOTONLY {

    // -- Run information update.

    cout << "Step: " << src -> step << "  Time: " << src -> time << endl;

    // -- Track particles.

    if (add) {
      Point          P, *I;
      FluidParticle* F;

      for (ListIterator<Point*> t (initial); t.more(); t.next()) {
	I   = t.current();
	P.x = I -> x;
	P.y = I -> y;
	P.z = I -> z;
	F   = new FluidParticle (src, FluidParticle::IDMax() + 1, P);
	if ((F -> inMesh())) particle.add (F);
      }
    }

    for (p.reset(); p.more(); p.next()) p.current() -> integrate ();
  }

  // -- CFL, energy, divergence information.

  if (cflstep && !(src -> step % cflstep)) {
    if (Geometry::nDim() == 3) modalEnergy();
    ROOTONLY {
      estimateCFL ();
      divergence  (work);
    }
  }

  // -- Periodic dumps and global information.
  
  const integer periodic = !(src->step %  (integer) Femlib::value("IO_HIS")) ||
                           !(src->step %  (integer) Femlib::value("IO_FLD"));
  const integer final    =   src->step == (integer) Femlib::value("N_STEP");
  const integer state    = periodic || final;

  if (state) {

    ROOTONLY {

      // -- Output particle locations.

      Point          P;
      FluidParticle* F;

      for (p.reset(); p.more(); p.next()) {
	F = p.current();
	if (F -> inMesh()) {
	  P = F -> location();
	  par_strm
	    << setw(10) << F -> ID()
	    << setw(15) << src -> time
	    << setw(15) << P.x
	    << setw(15) << P.y
	    << setw(15) << P.z
	    << endl;
	}
      }
    }

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

      ROOTONLY {
	his_strm << setw(4) << H->ID() << " " << setw(14) << src->time << " ";
	for (j = 0; j < NF; j++) his_strm << setw(15) << tmp[j];
	his_strm << endl;
      }
    }
     
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
  const integer    N     = Geometry::nModeProc();
  const integer    base  = Geometry::baseMode();
  const integer    nProc = Geometry::nProc();
  register integer i, m;
  vector<real>     ek (N);

  for (m = 0; m < N; m++) {
    ek[m] = 0.0;
    for (i = 0; i < DIM; i++) ek[m] += src -> u[i] -> mode_L2 (m);
  }

  if (nProc > 1) {

    ROOTONLY {
      for (m = 0; m < N; m++)
	mdl_strm << setw(10) << src -> time 
		 << setw( 5) << m 
		 << setw(15) << ek[m]
		 << endl;

      for (i = 1; i < nProc; i++) {
	Femlib::recv (ek(), N, i);
	for (m = 0; m < N; m++)
	  mdl_strm << setw(10) << src -> time 
		   << setw( 5) << m + i * N
		   << setw(15) << ek[m]
		   << endl;
      }
      
      mdl_strm.flush();

    } else
      Femlib::send (ek(), N, 0);

  } else 
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
  if (src -> nField() > 3) CFL_dt = max (CFL_dt, src -> u[2] -> CFL (2));

  dt_max  = SAFETY * CFL_max / CFL_dt;
  percent = (int) (100.0 * dt / dt_max);

  cout << "-- CFL: "     << CFL_dt * dt;
  cout << ", dt (max): " << dt_max;
  cout << ", dt (set): " << dt;
  cout << " ("           << percent << "%)" << endl;
}
