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

using namespace std;

#include <femdef.h>
#include <Utility.h>
#include <Femlib.h>

static char prog[] = "calc";
static void getargs (int, char**, istream*&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char     buf[StrMax];
  istream* input;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, input);

  while (input -> getline (buf, FILENAME_MAX))
    if (strstr (buf, "="))
      Femlib::value (buf);
    else
      cout << setprecision(17) << Femlib::value (buf) << endl;
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     istream*& input)
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
    input = new ifstream (*argv);
    if (input -> bad()) {
      cerr << usage;
      sprintf (buf, "unable to open file: %s", *argv);
      message (prog, buf, ERROR);
    }
  } else input = &cin;
}

