//////////////////////////////////////////////////////////////////////////////
// drive.C: compute solution to elliptic problem, optionally compare to
// exact solution.
//
// Copyright (C) 1994, 1999  Hugh Blackburn.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// USAGE:
// -----
// elliptic [options] session
//   options:
//   -h       ... print this message
//   -i       ... use iterative solver
//   -v[v...] ... increase verbosity level
//
//
// Author
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

#include <Sem.h>
#include <new.h>

static char prog[] = "elliptic";
static void memExhaust () { message ("new", "free store exhausted", ERROR); }
static void getargs    (int, char**, char*&);
static void getoptions (FEML*, char*&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, BoundarySys*&, Domain*&);

void Helmholtz (Domain*, const char*);


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
  
  char             *session, *forcing = 0, *exact = 0;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  BoundarySys*     bsys;
  Domain*          domain;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, session);

  preprocess (session, file, mesh, elmt, bman, bsys, domain);

  getoptions (file, forcing, exact);

  domain -> restart();

  Helmholtz (domain, forcing);

  ROOTONLY if (exact) domain -> u[0] -> errors (mesh, exact);

  domain -> dump();

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
  const char usage[] =
    "Usage: %s [options] session\n"
    "  options:\n"
    "  -h       ... print this message\n"
    "  -i       ... use iterative solver\n"
    "  -v[v...] ... increase verbosity level\n";
  char buf[StrMax];
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      ROOTONLY {
	sprintf (buf, usage, prog);
	cout << buf;
      }
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      Femlib::value ("ITERATIVE", (integer) 1);
      break;
    case 'v':
      do
	Femlib::value ("VERBOSE", (integer) Femlib::value ("VERBOSE") + 1);
      while (*++argv[0] == 'v');
      break;
    default:
      ROOTONLY {
	sprintf (buf, usage, prog);
	cout << buf;
      }
      exit (EXIT_FAILURE);
      break;
    }
  
  if (argc != 1) message (routine, "no session definition file", ERROR);

  session = *argv;
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


static void getoptions (FEML*  feml ,
			char*& force,
			char*& exact)
// ---------------------------------------------------------------------------
// Try to load forcing function string and exact solution string from USER
// section of FEML file.  The section is not required to be present.
// 
// Expect something in the form:
// <USER>
// forcing 0
// exact   sin(TWOPI*x)*sinh(TWOPI*y)/sinh(TWOPI)
// </USER>
//
// Either or both of the two strings may be absent.
// ---------------------------------------------------------------------------
{
  char routine[] = "options";
  char s[StrMax];

  if (feml -> seek ("USER")) {
    feml -> stream().ignore (StrMax, '\n');

    while (feml -> stream() >> s) {
      if (strcmp (s, "</USER>") == 0) break;

      upperCase (s);
      if (strcmp (s, "FORCING") == 0)
	feml -> stream() >> (force = new char [StrMax]);
      else if (strcmp (s, "EXACT") == 0)
	feml -> stream() >> (exact = new char [StrMax]);
    }

    if (strcmp (s, "</USER>") != 0)
      message (routine, "couldn't sucessfully close <USER> section", ERROR);
  }
}
