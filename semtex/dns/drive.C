//////////////////////////////////////////////////////////////////////////////
// drive.C
//
// SYNOPSIS:
// --------
// Control spectral element unsteady incompressible flow solver.
//
// USAGE:
// -----
// ns [options] session
//   options:
//   -h       ... print usage prompt
//   -i[i]    ... use iterative solver for viscous [and pressure] steps
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
// hmb@dbce.csiro.au
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <NS.h>
#include <new.h>

static char prog[] = "ns";
static void memExhaust () { message ("new", "free store exhausted", ERROR); }
static void getargs    (int, char**, char*&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);
#ifndef __DECCXX
  ios::sync_with_stdio();
#endif

  Geometry::CoordSys system;
  char      *session, fields[StrMax];
  int       np, nz, nel;
  FEML*     F;
  Mesh*     M;
  BCmgr*    B;
  Domain*   D;
  Analyser* A;
  
  Femlib::prep ();
  getargs      (argc, argv, session);
  
  F = new FEML (session);

  nz = (int) Femlib::value ("N_Z");
  
  M = new Mesh     (*F);
  B = new BCmgr    (*F);

  nel    = M -> nEl();  
  np     =  (int) Femlib::value ("N_POLY");
  nz     =  (int) Femlib::value ("N_Z"   );
  system = ((int) Femlib::value ("CYLINDRICAL") ) ?
                     Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, system);
  if   (nz > 1) strcpy (fields, "uvwp");
  else          strcpy (fields, "uvp");

  D = new Domain   (*F, *M, *B, fields, session);
  A = new Analyser (*D, *F);

  D -> initialize();

  NavierStokes (D, A);

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
  char routine[] = "getargs";
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
  
  if   (argc != 1) message (routine, "no session definition file", ERROR);
  else             session = *argv;
}
