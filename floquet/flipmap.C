///////////////////////////////////////////////////////////////////////////////
// flipmap.C: generate a list of index pairs for a symmetry-related transform.
//
// Copyright (C) 2002,2004 Hugh Blackburn.
//
// USAGE
// flipmap [options] [file]
// options:
//   -h       ... print this message
//   -x || -y ... reflection symmetry to be used to generate map
//   -t <num> ... set positional tolerance to num
//
// FILES
// Input file is a semtex/prism ASCII mesh file, generated e.g. by meshpr.
// Output is ASCII, with three header lines followed by a list of index pairs:
//   8 8 1 108 NR NS NZ NEL  # -- matches input file header
//   x                       # -- reflection symmetry generator
//   202                     # -- NFLIP, no. of pairs to follow < NR*NS*NEL/2
//   11   4001
//   12   4000
//   13   3999
//   ... etc. (NFLIP lines)
//
// The mapping between the mesh points may not be unique; we take the
// first mapping found for each case (reflect positive->negative,
// *and* negative->positive). This way a single gather will do the
// reflection, and leave no holes.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "sem.h"

static char prog[] = "flipmap";

const real EPS = EPSSP;	// -- Default positional tolerance.

static void    getargs  (int, char**, char&, real&, ifstream&);
static integer header   (ifstream&);
static integer loadmesh (ifstream&, const integer, const char,
			 vector<real>&, vector<real>&);
static void    findmap  (const char, const integer, const real tol,
			 const vector<real>&, const vector<real>&,
			 const integer, vector<integer>&, vector<integer>&);
static void     printup  (const char, const integer,
			  const vector<integer>&, const vector<integer>&);


int main (int    argc,
          char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  ifstream        file;
  char            generator;
  vector<real>    x, y;
  vector<integer> orig, flip;
  integer         npts, nmap;
  real            tol = EPS;

  getargs (argc, argv, generator, tol, file);

  npts = header   (file);
  nmap = loadmesh (file, npts, generator, x, y);
  
  findmap (generator, npts, tol, x, y, nmap, orig, flip);
  printup (generator, nmap, orig, flip);

  return EXIT_SUCCESS;
}


static void getargs (int       argc,
		     char**    argv,
		     char&     gen ,
		     real&     tol ,
		     ifstream& file)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: flipmap [-h] -x || -y [-t <num>] [meshfile]\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 't':
      if (*++argv[0]) tol = atof (  *argv);
      else { --argc;  tol = atof (*++argv); }
      break;
    case 'x':
      gen = 'x';
      break;
    case 'y':
      gen = 'y';
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (!(gen == 'x' || gen == 'y')) {
    cerr << prog << ": must specify -x or -y" << endl;
    exit (EXIT_FAILURE);
  }

  if (argc == 1)
    file.open (argv[0], ios::in);
  else {
    cerr << usage;
    exit (EXIT_FAILURE);
  }

  if (!file) {
    cerr << prog << ": unable to open file" << endl;
    exit (EXIT_FAILURE);
  }
}


static integer header (ifstream& file)
// ---------------------------------------------------------------------------
// If file header indicates this is a mesh file, output header and
// return number of points. Else die.
// ---------------------------------------------------------------------------
{
  integer  nr, ns, nz, nel;
  char     buf[StrMax];
  
  file >> nr >> ns >> nz >> nel;
  file.getline (buf, StrMax);

  if (!strstr (buf, "NR NS NZ NEL")) {
    cerr << prog
	 << ": mesh header line should include NR NS NZ NEL: " 
	 << buf << endl;
    exit (EXIT_FAILURE);
  }
  
  cout << nr << " " << ns << " " << 1 << " " << nel << " NR NS NZ NEL" << endl;

  return nr * ns * nel;
}


static integer loadmesh (ifstream&     file,
			 const integer npts,
			 const char    gen ,
			 vector<real>& x   ,
			 vector<real>& y   )
// ---------------------------------------------------------------------------
// Load the vectors x and y. Return the number of points to be mapped
// by symmetry.
// ---------------------------------------------------------------------------
{
  integer i, nmap = 0;

  x.resize (npts);
  y.resize (npts);

  for (i = 0; file && i < npts; i++) file >> x[i] >> y[i];

  if (!file) {
    cerr << prog << ": premature end of input" << endl;
    exit (EXIT_FAILURE);
  }

  if (gen == 'x') {
    for (i = 0; i < npts; i++) if (x[i] > EPS) nmap++;
  } else {
    for (i = 0; i < npts; i++) if (y[i] > EPS) nmap++;
  }

  nmap *= 2;

  if (nmap == 0) {
    cerr << prog << ": no points to map" << endl;
    exit (EXIT_FAILURE);
  }

  if (nmap > npts) {
    cerr << prog << ": too many points to map: "
	 << nmap << " vs. " << npts << endl;
    exit (EXIT_FAILURE);
  }

  return nmap;
}


static void findmap (const char          gen ,
		     const integer       npts,
		     const real          tol ,
		     const vector<real>& x   ,
		     const vector<real>& y   ,
		     const integer       nmap,
		     vector<integer>&    orig,
		     vector<integer>&    flip)
// ---------------------------------------------------------------------------
// This is where the mapping gets constructed. Order npts*npts operation.
// Take the first available mapping index for each point.
// ---------------------------------------------------------------------------
{
  integer i, j, k = 0;
  bool    found;

  orig.resize (nmap);
  flip.resize (nmap);

  if (gen == 'x') {
    for (i = 0; i < npts; i++) {
      if (x[i] > tol || x[i] < -tol) {
	for (found = false, j = 0; !found && j < npts; j++) {
	  if (fabs(x[i] + x[j]) < tol && fabs(y[i] - y[j]) < tol) {
	    orig[k]   = i;
	    flip[k++] = j;
	    found     = true;
	  }
	}
	if (!found) {
	  cerr << "Warning: mesh point "
	       << x[i] << ",\t" << y[i] << "\tnot mirrored" << endl;
	}
      }
    }
  } else {
    for (i = 0; i < npts; i++) {
      if (y[i] > tol || y[i] < -tol) {
	for (found = false, j = 0; !found && j < npts; j++) {
	  if (fabs(y[i] + y[j]) < tol && fabs(x[i] - x[j]) < tol) {
	    orig[k]   = i;
	    flip[k++] = j;
	    found     = true;
	  }
	}
	if (!found) {
	  cerr << "Warning: mesh point "
	       << x[i] << ",\t" << y[i] << "\tnot mirrored" << endl;
	}
      } 
    }
  }

  if (k != nmap) {
    cerr << prog
	 << ": number of maps found, " << k
	 << ", not equal to number needed, " << nmap << endl;
    exit (EXIT_FAILURE);
  }
}


static void printup (const char             gen ,
		     const integer          nmap,
		     const vector<integer>& orig,
		     const vector<integer>& flip)
// ---------------------------------------------------------------------------
// Print the map on standard output.
// ---------------------------------------------------------------------------
{
  integer i;

  cout << gen  << endl;
  cout << nmap << endl;
  
  for (i = 0; i < nmap; i++)
    cout << orig[i] << '\t' << flip[i] << endl;
}

