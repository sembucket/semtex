//////////////////////////////////////////////////////////////////////////////
// drive.C
//
// SYNOPSIS:
// --------
// Compute solution to elliptic problem, (compare to exact solution).
//
// USAGE:
// -----
// elliptic [options] session
//   options:
//   -h        ... print usage prompt
//   -i        ... use iterative solver
//   -v[v...]  ... increase verbosity level
//
//
// Author
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

#include <Sem.h>
#include <new.h>

static char prog[] = "elliptic";
static void memExhaust () { message ("new", "free store exhausted", ERROR); }

extern void Helmholtz (Domain*, const char*);
static void getargs   (int, char**, char*&);
static void setup     (FEML&, char*&, char*&);


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
  
  char* session;
  char* forcing = 0;
  char* exact   = 0;
  char  fields[StrMax];

  Femlib::prep();

  getargs (argc, argv, session);
  strcpy  (fields, "u");

  FEML  F (session);

  setup (F, forcing, exact);

  Mesh   M (F);
  BCmgr  B (F);
  Domain* D = new Domain (F, M, B, fields, session);

  D -> initialize();

  Helmholtz (D, forcing);

  if (exact) D -> u[0] -> errors (M, exact);

  char     outname[StrMax];
  ofstream output (strcat (strcpy (outname, F.root()), ".fld"));
  D -> dump  (output);

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
  int  level;
  char usage[] =
    "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -i       ... use iterative solver\n"
    "  -v[v...] ... increase verbosity level\n";
 
  while (--argc  && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      Femlib::value ("ITERATIVE", 1);
      break;
    case 'v':
      do
	Femlib::value ("VERBOSE", (int) Femlib::value ("VERBOSE") + 1);
      while (*++argv[0] == 'v');
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


static void setup (FEML&  feml ,
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
  char s[StrMax], g[StrMax], err[StrMax];

  if (feml.seek ("USER")) {
    feml.stream().ignore (StrMax, '\n');

    while (feml.stream() >> s) {
      if (strcmp (s, "</USER>") == 0) break;
      upperCase (s);

      if (strcmp (s, "FORCING") == 0) {
	force = new char [StrMax];
	feml.stream() >> force;

      } else if (strcmp (s, "EXACT") == 0) {
	exact = new char [StrMax];
	feml.stream() >> exact;
	
      } else if (strcmp (s, "GEOMETRY") == 0) {
	feml.stream() >> g;
	upperCase (g);
	if (strcmp (g, "2D-CARTESIAN")) {
	  sprintf (err, "bad geometry in USER section: %s", g);
	  message (routine, err, ERROR);
	} else
	  Femlib::value ("GEOMETRY", CART2D);

      } else {
	sprintf (err, "undefined in USER section: %s", s);
	message (routine, err, ERROR);
      }
    }

    if (strcmp (s, "</USER>") != 0)
      message (routine, "couldn't sucessfully close <USER> section", ERROR);
  }
}
