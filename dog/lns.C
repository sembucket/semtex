//////////////////////////////////////////////////////////////////////////////
// lns.C: control spectral element DNS for incompressible flows.
// This version drives linear evolution of a single Fourier mode.
//
// Copyright (C) 1994,2004 Hugh Blackburn & John Elston
//
// USAGE:
// -----
// lns [options] session
//   options:
//   -h       ... print usage prompt
//   -i[i]    ... use iterative solver for viscous [and pressure] steps
//   -v[v...] ... increase verbosity level
//   -chk     ... checkpoint field dumps
//
// AUTHOR:
// ------
// Hugh Blackburn
// CSIRO
// P.O. Box 56
// Highett, Vic 3190
// Australia
// hugh.blackburn@csiro.au
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "stab.h"

static char prog[] = "lns";
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
  BCmgr*           bman;
  BoundarySys*     bsys;
  Domain*          domain;
  StabAnalyser*    analyst;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, session);

  preprocess (session, file, mesh, elmt, bman, bsys, domain);
 
  analyst = new StabAnalyser (domain, file);

  domain -> restart ();
  domain -> loadBase();
  domain -> report  ();
  
  integrate (domain, analyst);

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
  const char routine[] = "getargs";
  char       buf[StrMax];
  char       usage[]   = "Usage: %s [options] session-file\n"
    "[options]:\n"
    "-h        ... print this message\n"
    "-i[i]     ... use iterative solver for viscous [& pressure] steps\n"
    "-v[v...]  ... increase verbosity level\n"
    "-chk      ... checkpoint field dumps\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      do
	Femlib::ivalue ("ITERATIVE", Femlib::ivalue ("ITERATIVE") + 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE", Femlib::ivalue ("VERBOSE") + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv)) {
	Femlib::ivalue ("CHKPOINT",  static_cast<int>(1));
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
			BCmgr*&           bman   ,
			BoundarySys*&     bsys   ,
			Domain*&          domain )
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed in order of creation.
// ---------------------------------------------------------------------------
{
  const integer verbose = Femlib::ivalue ("VERBOSE");
  const real*   z;
  integer       i, np, nel, npert;

  // -- Set default additional tokens.

  Femlib::value ("BASE_PERIOD", 0.0);

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  np    = Femlib::ivalue ("N_POLY");
  nel   = mesh -> nEl();
  npert = file -> attribute ("FIELDS", "NUMBER") - 1;
  
  Geometry::set (nel, npert);

  VERBOSE cout << "done" << endl;

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
}
