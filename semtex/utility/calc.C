///////////////////////////////////////////////////////////////////////////////
// calc: a basic calculator using the function parser.
//
// Copyright (c) 1994 Hugh Blackburn
//
// Usage: calc [-h] [file]
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <femdef.h>
#include <Utility.h>
#include <Femlib.h>

static char prog[] = "calc";
static void getargs (int, char**);


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  char buf[StrMax];

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv);

  while (cin.getline (buf, FILENAME_MAX))
    if (strstr (buf, "="))
      Femlib::value (buf);
    else
      cout << setprecision(17) << Femlib::value (buf) << endl;
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int    argc,
		     char** argv)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char buf[StrMax];
  char usage[] = "Usage: %s [-h] [file]\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cerr << "-- Calculator operators and functions:" << endl;
      yy_help ();
      cerr << endl << "-- Preset internal variables:"  << endl;
      yy_show ();
      exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    ifstream* inputfile = new ifstream (*argv);
    if (inputfile -> good()) {
      cin = *inputfile;
    } else {
      sprintf (buf, usage, prog);
      cerr << buf;
      sprintf (buf, "unable to open file: %s", *argv);
      message (prog, buf, ERROR);
    }
  }
}

