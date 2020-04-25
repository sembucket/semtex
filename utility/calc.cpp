/*****************************************************************************
 * calc: a simple calculator based on the femlib function parser.
 *
 * Usage
 * -----
 * calc [-h] [file]
 *
 * @file utility/calc.cpp
 * @ingroup group_utility
 *****************************************************************************/
// Copyright (c) 1994 <--> $Date: 2020/01/06 04:35:44 $, Hugh Blackburn
// --
// This file is part of Semtex.
// 
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: calc.cpp,v 9.2 2020/01/06 04:35:44 hmb Exp $"; 

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include <cfemdef.h>
#include <utility.h>
#include <femlib.h>

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
      cerr << "-- Preset internal variables:"  << endl;
      yy_show ();
      cerr << endl;
      cerr << "-- Calculator operators, functions and procedures:" << endl;
      yy_help ();
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
    if (input -> fail()) {
      cerr << usage;
      sprintf (buf, "unable to open file: %s", *argv);
      message (prog, buf, ERROR);
    }
  } else input = &cin;
}

