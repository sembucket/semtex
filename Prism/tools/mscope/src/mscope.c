/*
 *  Mscope (c) 1994-99 R. D. Henderson and H. M. Blackburn
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *                      *    *    *    *    *    *
 *
 *  Welcome to Mscope!  This is a program designed to examine, debug, refine,
 *  and otherwise pick apart a spectral element mesh.  It uses a simple
 *  graphical display to show mesh data and a command-line interface to 
 *  accept commands from the user. 
 *
 *  Mscope understands both conforming and nonconforming meshes, performs
 *  various automatic tests on the quality of the mesh and the numbering
 *  scheme, and can even indicate when you might have a bug in the mesh.  
 *  It also supports mesh generation by allowing the user to interactively
 *  add elements to a mesh.  This can be automated by directing mscope to
 *  read a script of mesh generation commands.
 *
 *  Mscope's primary purpose is to display element, patch, and node numbering 
 *  schemes graphically so that you can get a quick look at these things
 *  without staring at a blur of numbers on the screen.
 *
 *  Mscope uses a variant of Prism's spectral element library called "cubit".
 *  This supports the same basic data types as speclib but adds the function-
 *  ality needed for adaptive mesh refinement.
 * 
 *  Please send comments or suggestions to:
 *
 *  Ron Henderson
 *  Aeronautics and Applied Mathematics
 *  California Institute of Technology
 *  ron@galcit.caltech.edu
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <strings.h>
#include <stdlib.h>

#include "mscope.h"
#include "pl/pl.h"

static char *prog   = "mscope";
static char *usage  = "[options] file";
static char *help   = 
"options:\n"
"-b         ... run in batch mode w/out interactive graphics\n"
"-r file    ... specify a mesh file to read\n"
"-s script  ... process a script before switching to command mode\n"
"-d device  ... specify the graphics device (default = x11).  See mscope's\n"
"               internal help command for further options.\n";

static char *legal  =
"Mscope comes with ABSOLUTELY NO WARRANTY; for details type 'warranty'. This\n"
"is free software, and you are welcome to redistribute it under certain     \n"
"conditions; type 'copyright' for details.\n\n";

static char *rcsid  = "$Revision$";
static char *prompt = "Mscope> ";
static char *dev    = NULL;
static char *script = NULL;
static char *mesh   = NULL;

/* This is the mesh and field info shared by all application routines. */

Domain Geometry;
FILE*  mscope_command_stream = stdin;

/* ------------------------------------------------------------------------- */

int main (int argc, char *argv[])
{
  char buf[BUFSIZ];
  char *op;

  /* Initialize the libraries */

  pl_init();
  cubit_init();

  /* Initialize options and flags */

  option_set ("autoRedraw",  1);
  option_set ("autoScale",   1);
  option_set ("autoConnect", 1);
  option_set ("autoProfile", 1);
  option_set ("nekton",      M_NEKTON);
  option_set ("prism",       M_NEKTON);
  option_set ("semtex",      M_SEMTEX);

  Family_disable();

  /* Initialize options from environment variables */

  if ((op=getenv("MSCOPE_FORMAT")))
    option_set("format", option(op));

  /* Install library interfaces */

  user_refine = mscope_refine;
  user_bc     = mscope_bc;
  user_prune  = NULL;
  user_perm   = NULL;

  /* Parse command line arguments */

  parse_args (argc, argv);

  if (!option("batch")) {
    sscanf (rcsid, "%*s %s", buf);
    printf ("Mscope v%s\n" 
	    "Copyright (c) 1994-1999 "
	    "R. D. Henderson and H. M. Blackburn\n", buf);
    fputs  (legal, stdout);
    
    sprintf (buf, "dev %s\n", dev ? dev : "x11 -t Mscope");
    DoParse (buf);
  } else {
    sprintf (buf, "dev %s\n", dev ? dev : "nodevice");
    DoParse (buf);
  }
  
  /* ----------  Command Loop ---------- */
  
  if (mesh) {
    sprintf (buf, "meshin %s", mesh);
    DoParse (buf);
  }
  
  if (script) {
    sprintf (buf, "input %s", script);
    DoParse (buf);
  }

  if (!option("batch")) {
    do {
      fputs (prompt, stdout);
      fgets (buf, BUFSIZ, stdin);
    } 
    while (DoParse(buf) != EOF);
  }
  
  /* ----------  Command Loop ---------- */
  
  /* Make sure the active device is flushed before exiting */
  sprintf (buf, "dev nodevice");
  DoParse (buf);

  cubit_exit();
  return 0;
}

/* Parse command line arguments */

void parse_args (int argc, char *argv[])
{
  int n;

  for (n = 1; n < argc; n++) {
    if (argv[n][0] == '-') {
      const char c = argv[n][1];
      switch (c) {
      case 'b':
	option_set("batch", 1);
	break;
      case 'd':
	dev = strdup(argv[++n]);
	break;
      case 'h':
	fprintf (stderr, "usage: %s %s\n%s", prog, usage, help);
	exit (0);
      case 's':
	script = argv[++n];
	break;
      case 'r':
	mesh = argv[++n];
	break;
      default:
	fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }
      
      /* Try the various default suffixes */

    } 

    else if (strstr(argv[n],".rea"))
      mesh = strdup(argv[n]);
    else if (strstr(argv[n],".feml"))
      mesh = strdup(argv[n]);
    else if (strstr(argv[n],".ms"))
      script = strdup(argv[n]);
    else
      mesh = strdup(argv[n]);
  }
}

