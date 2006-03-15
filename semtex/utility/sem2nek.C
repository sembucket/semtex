///////////////////////////////////////////////////////////////////////////////
// sem2nek.C: read a semtex input file, print a NEKTON-style input file.
//
// Copyright (c) 1997 <--> $Date$,  Hugh Blackburn
//
// Usage:
// -----
// sem2nek [options] session
//   options:
//   -h       ... print usage prompt
//
// Files:
// -----
// The GROUPS BCS NODES ELEMENTS blocks of a semtex input file are required,
// the others need not be set.  Output concentrates on mesh & BC descriptions,
// other sections may need to be edited by hand.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <cstdlib>
#include <iomanip>

#include <cfemdef.h>
#include <femlib.h>
#include <mesh.h>
#include <utility.h>

static char prog[] = "sem2nek";

static void getArgs   (int, char**, char*&);
static void printHead ();
static void printTail ();


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  ifstream  file;
  char*     session;

  Femlib::initialize (&argc, &argv);
  getArgs (argc, argv, session);

  FEML feml (session);
  Mesh M (&feml);

  printHead ();
  M.printNek();
  printTail ();

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getArgs (int    argc   ,
		     char** argv   ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char routine[] = "getArgs";
  char buf[StrMax];
  char usage[]   =
    "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h        ... print this message\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      message ("", buf, REMARK);
      exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (buf, usage, prog);
      message ("", buf, REMARK);
      exit (EXIT_FAILURE);
      break;
    }
  
  if   (argc != 1) message (routine, "no session definition file", ERROR);
  else             session = *argv;
}


static void printHead ()
// ---------------------------------------------------------------------------
// Print the beginning of a NEKTON .rea file.  May need editing.
// ---------------------------------------------------------------------------
{
  int  np = (int) Femlib::value ("N_P");
  int  nz = (int) Femlib::value ("N_Z");
  real dt =       Femlib::value ("D_T");

  cout << "****** PARAMETERS *****"                                  << endl;
  cout << "6 PRISM VERSION"                                          << endl;
  cout << "2 DIMENSIONAL RUN"                                        << endl;
  cout << "8 PARAMETERS FOLLOW"                                      << endl;
  cout << "100.0" "\tRe"     "\t# Reynolds Number"                   << endl;
  cout << "1"     "\tEQTYPE" "\t# Skew-symmetric form"               << endl;
  cout << np <<   "\tNORDER" "\t# Polynomial order + 1"              << endl;
  cout << dt <<   "\tDT"     "\t# Time step"                         << endl;
  cout << "20*DT" "\tTIME"   "\t# Integration time"                  << endl;
  cout << "0.02"  "\tIO_FLD" "\t# I/O flush & checkpoint"            << endl;
  cout << "10"    "\tIO_HIS" "\t# Measurements of lift, drag, etc."  << endl;
  cout << "50"    "\tIO_CFL" "\t# CFL estimate"                      << endl;
  cout << "0"     "\tLines of passive scalar data follows"           << endl;
  cout << "0"     "\tLOGICAL SWITCHES FOLLOW"                        << endl;
}


static void printTail ()
// ---------------------------------------------------------------------------
// Print the end of a NEKTON .rea file.  May need editing.
// ---------------------------------------------------------------------------
{
  cout << "***** NO THERMAL BOUNDARY CONDITIONS *****"                << endl;
  cout << "1     INITIAL CONDITIONS *****"                            << endl;
  cout << "Default"                                                   << endl;
  cout << "***** DRIVE FORCE DATA ***** PRESSURE GRAD, FLOW, Q"       << endl;
  cout << "0                 Lines of Drive force data follow"        << endl;
  cout << "***** Variable Property Data *** Overrides Parameter data" << endl;
  cout << "1 Lines follow."                                           << endl;
  cout << "0 PACKETS OF DATA FOLLOW"                                  << endl;
  cout << "***** HISTORY AND INTEGRAL DATA *****"                     << endl;
  cout << "0   POINTS.  Hcode, I,J,H,IEL"                             << endl;
  cout << "***** OUTPUT FIELD SPECIFICATION *****"                    << endl;
  cout << "6 SPECIFICATIONS FOLLOW"                                   << endl;
  cout << "F      COORDINATES"                                        << endl;
  cout << "T      VELOCITY"                                           << endl;
  cout << "T      PRESSURE"                                           << endl;
  cout << "F      TEMPERATURE"                                        << endl;
  cout << "F      TEMPERATURE GRADIENT"                               << endl;
  cout << "0      PASSIVE SCALARS OUTPUTS"                            << endl;
  cout << "***** OBJECT SPECIFICATION *****"                          << endl;
  cout << "0 Surface Objects"                                         << endl;
  cout << "0 Volume  Objects"                                         << endl;
  cout << "0 Edge    Objects"                                         << endl;
  cout << "0 Point   Objects"                                         << endl;
  cout << "***** CURVED SIDE FORTRAN FUNCTIONS *****"                 << endl;
  cout << "0 Lines of Fortran"                                        << endl;
}

