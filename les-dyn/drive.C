//////////////////////////////////////////////////////////////////////////////
// drive.C: control spectral element LES for incompressible flows.
//
// This version for dynamic eddy-viscosity based SGSS, with Smagorinsky
// as the underlying model.
//
// Copyright (c) 1999--2000 Hugh Blackburn
//
// USAGE:
// -----
// les-dyn [options] session
//   options:
//   -h       ... print usage prompt
//   -v[v...] ... increase verbosity level
//
// AUTHOR:
// ------
// Hugh Blackburn
// CSIRO Division of Building, Construction and Engineering
// P.O. Box 56
// Highett, Vic 3190
// Australia
// hugh.blackburn@dbce.csiro.au
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "les.h"
#include <new.h>

static char prog[] = "les-dyn";
static void memExhaust () { message ("new", "free store exhausted", ERROR); }
static void getargs      (int, char**, char*&);
static void preprocess   (const char*, FEML*&, Mesh*&, vector<Element*>&,
			  BCmgr*&, BoundarySys*&, Domain*&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);
#if !defined(__alpha)
  ios::sync_with_stdio();
#endif

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
  
  integrate (domain, new LESAnalyser (domain, file));

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
	Femlib::value ("ITERATIVE", (integer) Femlib::value ("ITERATIVE") + 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::value ("VERBOSE",   (integer) Femlib::value ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv)) {
	Femlib::value ("CHKPOINT",  (integer) 1);
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
  const integer      verbose = (integer) Femlib::value ("VERBOSE");
  Geometry::CoordSys space;
  const real*        z;
  integer            i, np, nz, nel;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   = mesh -> nEl();
  np    =  (integer) Femlib::value ("N_POLY");
  nz    =  (integer) Femlib::value ("N_Z");
  space = ((integer) Femlib::value ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;

  if (nz < 2) message (prog, "3D only, N_Z > 1 required", ERROR);
  
  Geometry::set (np, nz, nel, space);

  VERBOSE cout << "done" << endl;

  // -- Build all the elements.

  VERBOSE cout << "Building elements ... ";

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  elmt.setSize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, mesh, z, np);

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
