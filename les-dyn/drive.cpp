//////////////////////////////////////////////////////////////////////////////
// drive.C: control spectral element LES for incompressible flows.
//
// This version for dynamic eddy-viscosity based SGSS, with Smagorinsky
// as the underlying model.
//
// Copyright (c) 1999, 2001 Hugh Blackburn, Stefan Schmidt
//
// USAGE:
// -----
// les-dyn [options] session
//   options:
//   -h       ... print usage prompt
//   -v[v...] ... increase verbosity level
//
// REFS:
// ----
//
// @Article{blsc03,
// author = 	 {H. M. Blackburn and S. Schmidt},
// title = 	 {Spectral Element Filtering Techniques for Large Eddy
//                 Simulation with Dynamic Estimation},
// journal = 	 JCP,
// year = 	 2003,
// volume =       186,
// number =       2,
// pages =        {610--629}
// }
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

#include "les.h"

static char prog[] = "les-dyn";
static void getargs    (int, char**, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, BoundarySys*&, Domain*&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char*            session;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bmgr;
  BoundarySys*     bsys;
  Domain*          domain;
  
  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session);

  preprocess  (session, file, mesh, elmt, bmgr, bsys, domain);
  initFilters ();

  domain -> restart();

  ROOTONLY domain -> report();
  
  integrate (domain, new LESAnalyser(domain, file), new SumIntegrator(domain));

  Femlib::finalize();

  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  const char *routine = "getargs";
  char       buf[StrMax], *usage =
    "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -v[v...] ... increase verbosity level\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      do
	Femlib::value ("ITERATIVE", (int_t) Femlib::value ("ITERATIVE") + 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::value ("VERBOSE",   (int_t) Femlib::value ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv)) {
	Femlib::value ("CHKPOINT",  (int_t) 1);
      } else {
	fprintf (stdout, usage, prog);
	exit (EXIT_FAILURE);	  
      }
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if   (argc != 1) message (routine, "no session definition file", ERROR);
  else             session = *argv;
}


static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bmgr   ,
			BoundarySys*&     bsys   ,
			Domain*&          domain )
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed in order of creation.
// ---------------------------------------------------------------------------
{
  const int_t        verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  const real_t*      z;
  int_t              i, np, nz, nel;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   = mesh -> nEl();
  np    =  Femlib::ivalue ("N_P");
  nz    =  Femlib::ivalue ("N_Z");
  space = (Femlib::ivalue ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;

  if (nz < 2) message (prog, "3D only, N_Z > 1 required", ERROR);
  
  Geometry::set (np, nz, nel, space);

  VERBOSE cout << "done" << endl;

  // -- Build all the elements.

  VERBOSE cout << "Building elements ... ";

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

  VERBOSE cout << "done" << endl;

  // -- Build all the boundary condition applicators.

  VERBOSE cout << "Building boundary condition manager ..." << endl;

  bmgr = new BCmgr (file, elmt);

  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.

  VERBOSE cout << "Building domain ..." << endl;

  domain = new Domain (file, elmt, bmgr);

  VERBOSE cout << "done" << endl;
}