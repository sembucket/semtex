//////////////////////////////////////////////////////////////////////////////
// drive.C
//
// SYNOPSIS:
// --------
// Compute solution to elliptic problem, optionally compare to exact solution.
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

extern void Helmholtz (Domain*, const char*);
static void getargs   (int, char**, char*&);
static void setup     (FEML*, char*&, char*&);


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
  
  char               *session, *forcing = 0, *exact = 0, field[2];
  integer            np, nz, nel;
  Geometry::CoordSys space;
  FEML*              F;
  Mesh*              M;
  BCmgr*             B;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, session);
  strcpy  (field, "c");

  F = new FEML  (session);
  M = new Mesh  (*F);
  B = new BCmgr (*F);

  nel   = M -> nEl();
  np    =  (integer) Femlib::value ("N_POLY");
  nz    =  (integer) Femlib::value ("N_Z");
  space = ((integer) Femlib::value ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);

  Domain* D = new Domain (*F, *M, *B, field, session);

  D -> initialize();

  setup (F, forcing, exact);

  Helmholtz (D, forcing);

  ROOTONLY if (exact) D -> u[0] -> errors (*M, exact);

  D -> dump();

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


static void setup (FEML*  feml ,
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
  char routine[] = "setup";
  char s[StrMax];

  if (feml -> seek ("USER")) {
    feml -> stream().ignore (StrMax, '\n');

    while (feml -> stream() >> s) {
      if (strcmp (s, "</USER>") == 0) break;
      upperCase (s);

      if (strcmp (s, "FORCING") == 0) {
	force = new char [StrMax];
	feml -> stream() >> force;

      } else if (strcmp (s, "EXACT") == 0) {
	exact = new char [StrMax];
	feml -> stream() >> exact;
	
      }
    }

    if (strcmp (s, "</USER>") != 0)
      message (routine, "couldn't sucessfully close <USER> section", ERROR);
  }
}
