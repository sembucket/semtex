//////////////////////////////////////////////////////////////////////////////
// drive.C:
//
// SYNOPSIS:
// --------
// Control spectral element aeroelastic flow solver.
//
// USAGE:
// -----
// aero [options] session
//   options:
//   -h       ... print usage prompt
//   -i[i]    ... use iterative solver for viscous [and pressure] steps
//   -v[v...] ... increase verbosity level
//   -chk     ... checkpoint field dumps
//
// FILES:
// -----
// Body motion parameters are set up in a file named session.bdy.
// See comments in body.C for file structure.
//
// AUTHOR:
// ------
// Hugh Blackburn
// CSIRO Division of Building, Construction and Engineering
// P.O. Box 56
// Highett, Vic 3190
// Australia
// hmb@dbce.csiro.au
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <aero.h>
#include <new.h>


static char prog[]  = "aero";
static void memExhaust () { message ("new", "free store exhausted", ERROR); }
static void getargs (int, char**, char*&);

void NavierStokes (Domain*, Body*, Analyser*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);
#ifndef __DECCXX
  ios::sync_with_stdio ();
#endif

  char      *session, fields[StrMax];
  int       nz;
  FEML*     F;
  Mesh*     M;
  BCmgr*    B;
  Domain*   D;
  Analyser* A;
  Body*     BD;

  Femlib::prep ();
  getargs      (argc, argv, session);

  cout << prog << ": aeroelastic Navier--Stokes solver"  << endl;
  cout << "      (c) Hugh Blackburn 1995--97."   << endl << endl;
  
  F = new FEML (session);

  nz = (int) Femlib::value ("N_Z");
  if (nz > 1) {
    if (nz & 1) {
      sprintf (fields, "N_Z (%1d) must be even", nz);
      message (prog, fields, ERROR);
    }
    strcpy (fields, "uvwp");
  } else
    strcpy (fields, "uvp");
  
  M  = new Mesh    (*F);
  B  = new BCmgr   (*F);
  D  = new Domain  (*F, *M, *B, fields, session);
  BD = new Body    (session);

  D  -> initialize();
  BD -> force (*D);

  A = new Analyser (*D, *BD);

  NavierStokes (D, BD, A);

#if 0
  ifstream*  input = new ifstream;
  char*      session;
  char       s[StrMax];

  // -- Initialization section.

//  FamilyMgr::active = 0;

  cout << prog << ": aeroelastic Navier--Stokes solver"  << endl;
  cout << "      (c) Hugh Blackburn 1995, 1996." << endl << endl;

  Femlib::prep  ();
  getArgs       (argc, argv, session);
  input -> open (session);
  setUp         (*input);

  // -- Get mesh information, save BC and parameter information.

  Mesh*  M = preProcess (*input);
  input -> close ();
  
  // -- Set up domain with single field, 'u'.

  Domain*  D = new Domain (*M, session, Femlib::integer ("N_POLY"));
  D -> u[0] -> setName ('u');

  // -- Add remaining velocity fields.

  const int DIM = Femlib::integer ("N_VAR");
  SystemField*  newField;
  for (int i = 1; i < DIM; i++) {
    newField = new SystemField (*D -> u[0]);
    D -> addField (newField);
    D -> u[i] -> setName ('u' + i);
  }

  // -- And constraint field 'p'.

  SystemField* Pressure = new SystemField (*D -> u[0], 1);
  D -> addField (Pressure);
  Pressure -> setName ('p');  
  PBCmanager::build (*Pressure);
  Pressure -> connect (*M, Femlib::integer ("N_POLY"));

  // -- Startup.

  D -> restart ();

  // -- Seek body information, construct body.

  input -> open (strcat (strcpy (s, session), ".bdy"));
  Body*  B = new Body (*input);
  input -> close ();

  B -> force   (*D);

  Analyser*  A = new Analyser (*D, *B);

  // -- Solve.

  NavierStokes (D, B, A);
#endif

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
  char buf[StrMax], c;
  char usage[]   =
    "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h        ... print this message\n"
    "  -i[i]     ... use iterative solver for viscous [& pressure] steps\n"
    "  -v[v...]  ... increase verbosity level\n"
    "  -chk      ... checkpoint field dumps\n";
 
  while (--argc  && **++argv == '-')
    switch (c = *++argv[0]) {
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
  
  if   (argc != 1) message (prog, "no session definition file", ERROR);
  else             session = *argv;  
}
