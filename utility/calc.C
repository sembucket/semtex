///////////////////////////////////////////////////////////////////////////////
// calc: a basic calculator using the function parser.
//
// Copyright (c) 1994--1999 Hugh Blackburn
//
// Usage: calc [-h] [file]
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include <femdef.h>
#include <Utility.h>
#include <Femlib.h>

static char prog[] = "calc";
static void getargs (integer, char**, ifstream&);


integer main (integer argc,
	      char**  argv)
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  char     buf[StrMax];
  ifstream file;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, file);

  while (file.getline (buf, FILENAME_MAX))
    if (strstr (buf, "="))
      Femlib::value (buf);
    else
      cout << setprecision(17) << Femlib::value (buf) << endl;
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (integer   argc,
		     char**    argv,
		     ifstream& file)
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
  
  if   (argc == 1) file.open   (*argv, ios::in);
  else             file.attach (0);

  if (!file) {
    sprintf (buf, usage, prog);
    cerr << buf;
    sprintf (buf, "unable to open file: %s", *argv);
    message (prog, buf, ERROR);
  }
}

