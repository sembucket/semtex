//////////////////////////////////////////////////////////////////////////////
// drive.C: control spectral element solver for unsteady
// incompressible flow with scalar transport (which adds a extra
// advection--diffusion equation for concentration 'c').
//
// Define token PRANDTL = (kinematic viscosity)/( scalar diffusivity).
// Default value 0.720.
//
// Optional Boussinesq buoyancy. The Boussinesq buoyancy term in the
// momentum equation is
//                      - BETA_T * (T - T_REF) g.
// Tokens used for buoyancy:
//
// GRAVITY:       magnitude of gravity vector.
// g_1, g_2, g_3: direction cosines of the gravity vector (normalised to 1).
// FFC:           (optional) uniform forcing term for scalar equation.
// T_REF:         reference temperature.
// BETA_T:        coefficient of (thermal) expansion = T_REF^(-1) for a gas.

//
// NB: for cylindrical coordinates, only axial gravity is currently
// implemented.

// Copyright (C) 1997 <--> $Date$, Hugh Blackburn
//
// USAGE:
// -----
// scat [options] session
//   options:
//   -h       ... print usage prompt
//   -f       ... freeze velocity field (scalar advection/diffusion only)
//   -i[i]    ... use iterative solver for viscous [and pressure] steps
//   -t[t]    ... select time-varying BCs for mode 0 [or all modes]
//   -v[v...] ... increase verbosity level
//   -chk     ... turn off checkpoint field dumps [default: selected]
//
// AUTHOR:
// ------
// Hugh Blackburn
// Department of Mechanical & Aerospace Engineering
// Monash University
// Vic 3800
// Australia
// hugh.blackburn@eng.monash.edu.au
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <scat.h>

static char prog[] = "scat";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, BoundarySys*&, Domain*&);

void NavierStokes  (Domain*, BCmgr*, ScatAnalyser*);
void AdvectDiffuse (Domain*, BCmgr*, ScatAnalyser*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
#ifdef _GNU_SOURCE
  feenableexcept (FE_OVERFLOW);    // -- Force SIG8 crash on FP overflow.
#endif

  char*            session;
  bool             freeze = false;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  BoundarySys*     bsys;
  Domain*          domain;
  ScatAnalyser*    analyst;
  
  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);

  preprocess (session, file, mesh, elmt, bman, bsys, domain);

  analyst = new ScatAnalyser (domain, file);

  domain -> restart();

  ROOTONLY domain -> report();

  if (freeze) AdvectDiffuse (domain, bman, analyst);
  else        NavierStokes  (domain, bman, analyst);

  Femlib::finalize();

  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     bool&  freeze ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char buf[StrMax];
  char usage[]   =
    "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -f       ... freeze velocity field (scalar advection/diffusion only)\n"
    "  -i       ... use iterative solver for viscous step\n"
    "  -v[v...] ... increase verbosity level\n"
    "  -chk     ... turn off checkpoint field dumps [default: selected]\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'f':
      freeze = true;
      break;
    case 'i':
      do
	Femlib::ivalue ("ITERATIVE", 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE",   Femlib::ivalue ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv))     Femlib::ivalue ("CHKPOINT",    0);
      else { fprintf (stdout, usage, prog); exit (EXIT_FAILURE); }
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }
  
  if   (argc != 1) message (prog, "no session definition file", ERROR);
  else             session = *argv;

  Femlib::value ("DTBDX", 0.0);
}


static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			BoundarySys*&     bsys   ,
			Domain*&          domain )
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed in order of creation.
// ---------------------------------------------------------------------------
{
  const char routine[] = "preprocess";
  const int_t        verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  const real_t*      z;
  int_t              i, np, nz, nel, procid, seed;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   =  mesh -> nEl();
  np    =  Femlib::ivalue ("N_P");
  nz    =  Femlib::ivalue ("N_Z");
  space = (Femlib::ivalue ("CYLINDRICAL")) ?
                     Geometry::Cylindrical : Geometry::Cartesian;
   
  Geometry::set (np, nz, nel, space);

  VERBOSE cout << "done" << endl;

  // -- If token RANSEED > 0 then initialize the random number
  //    generator based on wall clock time and process ID (i.e. a "truly"
  //    pseudo-random number).  NB: it is important to have done this
  //    before any other possible call to random number routines.

  if (Femlib::ivalue("RANSEED") > 0) {
    procid = Geometry::procID();
    seed   = -abs((procid + 1) * (char) time(NULL));
  } else seed = -1;
  Veclib::ranInit (seed);

  // -- Build all the elements.

  VERBOSE cout << "Building elements ... ";

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

  VERBOSE cout << "done" << endl;

  // -- Build all the boundary condition applicators.

  VERBOSE cout << "Building boundary condition manager ..." << endl;

  bman = new BCmgr (file, elmt);

  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.

  VERBOSE cout << "Building domain ..." << endl;

  domain = new Domain (file, elmt, bman);

  VERBOSE cout << "done" << endl;

  // -- Sanity checks on installed tokens.

  if (Femlib::ivalue ("SVV_MN") > Geometry::nP())
    message (routine, "SVV_MN exceeds N_P", ERROR);
  if (Femlib::ivalue ("SVV_MZ") > Geometry::nMode())
    message (routine, "SVV_MZ exceeds N_Z/2", ERROR);
}
