//////////////////////////////////////////////////////////////////////////////
// drive.C
//
// Copyright (C) 1999 Hugh Blackburn
//
// SYNOPSIS:
// --------
// Control spectral element unsteady incompressible flow solver for
// simulation of advective/diffusive transport in coated capillary tube.
// Cartesian coordinates only.
//
// USAGE:
// -----
// chroma [options] mobile stationary
//   options:
//   -h       ... print usage prompt
//   -i       ... use iterative solver for diffusion step.
//   -v[v...] ... increase verbosity level
//   -chk     ... checkpoint field dumps
//
// AUTHOR:
// ------
// Hugh Blackburn
// CSIRO Division of Building, Construction and Engineering
// P.O. Box 56
// Highett, Vic 3190
// Australia
// hugh.blackburn@dbce.csiro.au
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include "Sem.h"
#include <new.h>

static char  prog[] = "chroma";
static void  memExhaust () { message ("new", "free store exhausted", ERROR); }
static void  getargs    (int, char**, char*&, char*&);
static void  preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			 BCmgr*&, BoundarySys*&, Domain*&);
static real* gasdata    (const char*, const int, const int);

void transport (Domain*,Domain*, Analyser*,Analyser*, const real*, MixPatch*);

#define SET1 Geometry::set (np, nz, nel1, space)
#define SET2 Geometry::set (np, nz, nel2, space)


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);
#if !defined(__DECCXX)
  ios::sync_with_stdio();
#endif

  char               *session1, *session2;
  FEML               *file1,    *file2;
  Mesh               *mesh1,    *mesh2;
  BCmgr              *bman1,    *bman2;
  BoundarySys        *bsys1,    *bsys2;
  Domain             *domain1,  *domain2;
  vector<Element*>   elmt1,     elmt2;
  Analyser           *IO1,      *IO2;
  real*              gas;
  MixPatch           *patch;

  integer            np, nz, nel1, nel2;
  Geometry::CoordSys space = Geometry::Cartesian;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, session1, session2);
  
  preprocess (session1, file1, mesh1, elmt1, bman1, bsys1, domain1);
  np   = Geometry::nP();
  nz   = Geometry::nZ();
  nel1 = Geometry::nElmt();

  IO1 = new Analyser (domain1, file1);

  preprocess (session2, file2, mesh2, elmt2, bman2, bsys2, domain2);
  nel2 = Geometry::nElmt();

  IO2 = new Analyser (domain2, file2);

  if (np != Geometry::nP()) message (prog, "polynomial order mismatch", ERROR);

  SET1; domain1 -> restart();
  SET2; domain2 -> restart();

  SET1; gas = gasdata (session1, np, nel1);

  patch = new MixPatch (domain1, domain2);

  transport (domain1, domain2, IO1, IO2, gas, patch);

  Femlib::finalize();

  return EXIT_SUCCESS;
}


static void getargs (int    argc    ,
		     char** argv    ,
		     char*& session1,
		     char*& session2)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char buf[StrMax];
  char usage[]   =
    "Usage: %s [options] mobile-phase-session stationary-phase-session\n"
    "  [options]:\n"
    "  -h        ... print this message\n"
    "  -i[i]     ... use iterative solver\n"
    "  -v[v...]  ... increase verbosity level\n"
    "  -chk      ... checkpoint field dumps\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      do
	Femlib::value ("ITERATIVE", (int) Femlib::value ("ITERATIVE") + 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::value ("VERBOSE",   (int) Femlib::value ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv)) {
	Femlib::value ("CHKPOINT", 1);
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
  
  if (argc != 2) message (prog, "no session definition files", ERROR);

  session1 = *  argv;
  session2 = *++argv;
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
  const integer      verbose = (integer) Femlib::value ("VERBOSE");
  Geometry::CoordSys space;
  const real*        z;
  integer            i, np, nz, nel;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (*file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   = mesh -> nEl();
  np    =  (integer) Femlib::value ("N_POLY");
  nz    =  (integer) Femlib::value ("N_Z");
  space = ((integer) Femlib::value ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;
  
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

  bman = new BCmgr (file, elmt);

  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.

  VERBOSE cout << "Building domain ..." << endl;

  domain = new Domain (file, elmt, bman);

  VERBOSE cout << "done" << endl;
}


static real* gasdata (const char* session,
		      const int   np     ,
		      const int   nel    )
// ---------------------------------------------------------------------------
// Load the (w-only) 2D gas phase velocity data from file
// "session.gas" This is a standard 2D SEM field file.  Check that it
// has the correct element order and number of elements matches input
// declaration nel.
// ---------------------------------------------------------------------------
{
  const char routine[] = "gasdata";
  char       buf[StrMax], err[StrMax];
  int        i, NP, NZ, NEL, swab;
  const int  nProc = Geometry::nProc();
  real*      w = new real [Geometry::planeSize()];

  Veclib::zero (Geometry::planeSize(), w, 1);

  ROOTONLY {
    strcat (strcpy (buf, session), ".gas");
    ifstream file  (buf);

    if (!file) {
      sprintf (err, "file with gas phase velocity data, %s, not found", buf);
      message (routine, err, ERROR);
    }

    file.getline(buf, StrMax);
    if (!strstr (buf, "Session")) {
      sprintf (err, "%s.gas doesn't look like a field file", session);
      message (routine, err, ERROR);
    }

    file.getline(buf, StrMax).getline(buf, StrMax);
    istrstream (buf, strlen(buf)) >> NP >> NP >> NZ >> NEL;
    
    if (NP != np) {
      sprintf (err, "size mismatch, np: %1d/%1d", np, NP);
      message (routine, err, ERROR);
    }

    if (NEL != nel) {
      sprintf (err, "size mismatch: nel-gas (%1d) != nel (%1d)", NEL, nel);
      message (routine, err, ERROR);
    }

    file.getline(buf, StrMax).getline(buf, StrMax).getline(buf, StrMax);
    file.getline(buf, StrMax).getline(buf, StrMax).getline(buf, StrMax);
    file.getline(buf, StrMax);

    Veclib::describeFormat (err);

    if (!strstr (buf, "binary"))
      message (routine, "input field file not in binary format", ERROR);
  
    if (!strstr (buf, "endian"))
      message (routine, "input field file in unknown binary format", WARNING);
    else {
      swab = ((strstr (buf, "big") && strstr (err, "little")) ||
	      (strstr (err, "big") && strstr (buf, "little")) );
    }

    file.read ((char*) w, Geometry::nPlane() * sizeof (real));
    if (swab) Veclib::brev (Geometry::nPlane(), w, 1, w, 1);

    cout << "-- Gas velocity field      : read from "
	 << session << ".gas" << endl;

    if (nProc > 1)
      for (i = 1; i < nProc; i++)
	Femlib::send (w, Geometry::planeSize(), i);
  } else			// -- This must be a parallel run.
    Femlib::recv (w, Geometry::planeSize(), 0);

  return w;
}
