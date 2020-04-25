/*****************************************************************************
 * data2df_template: a boilerplate utility file that could be a starting
 * point for further limited data processing.
 *
 * Usage
 * -----
 * data2df_template [options] [file]
 * options:
 * -h       ... print this message.
 *
 * Synopsis
 * --------
 * Read in a data file header and its data (no geometric information)
 * and output.  If file is not present, read from standard input.
 * Write to standard output. In this boilerplate file, no other action
 * is taken. The main intention is to show how to get data in and out.
 * Actions available in src/data2df.cpp (2D x Fourier) may be used.
 * 
 * @file utility/data2df_template.cpp
 * @ingroup group_utility
 *****************************************************************************/
// Copyright (c) 2008 <--> $Date: 2020/01/06 04:35:44 $, Hugh Blackburn

static char RCS[] = "$Id: data2df_template.cpp,v 9.2 2020/01/06 04:35:44 hmb Exp $";

#include <sem.h>
#include <data2df.h>

static char prog[] = "data2df_template";
static void getargs  (int, char**, istream*&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int_t            i;
  istream*         input;
  Header           h;
  vector<Data2DF*> u;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, input);

  *input >> h;
  cout   << h;

  u.resize (h.nFields());
  for (i = 0; i < h.nFields(); i++) {

    // - Create u[i].

    u[i] = new Data2DF (h.nr, h.nz, h.nel, h.flds[i]);

    // -- Read in u[i] and byte-swap if necessary.

    *input >> *u[i];
    if (h.swab()) u[i] -> reverse();

    // -- Do something to u[i].

    // -- Output u[i];

    cout << *u[i];
  }

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     istream*& input)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: data2df_template [options] [file]\n"
    "options:\n"
    "-h       ... print this message\n";
    
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> fail()) message (prog, "unable to open input file", ERROR);
  } else input = &cin;
}
