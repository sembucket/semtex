///////////////////////////////////////////////////////////////////////////////
// embed.C: from a file of ASCII data (optionally stdin), an integer
// data window length, and an integer skip/delay, construct the
// covariance matrix from the total number of lagged vectors that can
// be constructed, and find its eigenvalues and eigenvectors.
//
// Copyright (c) 2001, Hugh Blackburn.
//
// Usage:
// ------
// embed [-h] -n <num> -s <num> [-p <num>] [file]
//
// Options:
// --------
// -h      : print usage string.
// -n <num>: select length of data window/number of eigenpairs.
// -s <num>: select skip/delay in input data.
// -p <num>: project original input onto the p leading eigenvectors.
//
// Reference:
// ----------
// D.S. Broomhead and G. P. King (1986), Extracting qualitative
// dynamics from experimental data, Physica D 20, 217--236.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>

#include <Utility.h>
#include <Veclib.h>
#include <Stack.h>
#include <Array.h>
#include <Blas.h>
#include <Lapack.h>

static char prog[]  = "embed";
static char usage[] = "embed [-h] -n <num> -s <num> [-p <num>] [file]";

static void getargs  (int, char**, int&, int&, int&);
static void getdata  (istream&, vector<double>&, const int, const int);
static void covary   (vector<double>&, vector<double>&, const int);
static void eigensys (vector<double>&, vector<double>&, const int);
static void putsys   (const vector<double>&, const int);
static void project  (const vector<double>&, const int,
		      const vector<double>&, const int);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  vector<double> data, cov, esys;
  int            s = 0, n = 0, p = 0;

  getargs  (argc, argv, s, n, p);
  getdata  (cin, data, s, n);

  covary   (data, cov, n);
  eigensys (cov, esys, n);
  
  putsys   (esys, n);
  project  (data, p, esys, n);

  return EXIT_SUCCESS;
}


static void getargs (int       argc,
		     char**    argv,
		     int&      skip,
		     int&      wind,
		     int&      proj)
// ---------------------------------------------------------------------------
// Parse command line arguments. Skip and wind are mandatory.
// ---------------------------------------------------------------------------
{
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage << endl;
      exit (EXIT_SUCCESS);
      break;
    case 'n':
      if (*++argv[0])
	wind = atoi (*argv);
      else {
	wind = atoi (*++argv);
	argc--;
      }
      break;
    case 'p':
      if (*++argv[0])
	proj = atoi (*argv);
      else {
	proj = atoi (*++argv);
	argc--;
      }
      break;
    case 's':
      if (*++argv[0])
	skip = atoi (*argv);
      else {
	skip = atoi (*++argv);
	argc--;
      }
      break;
    default:
      cerr << usage << endl;
      exit (EXIT_FAILURE);
      break;
    }

  if (!skip || !wind) { cerr << usage << endl; exit (EXIT_FAILURE); }
  if (proj < 0 || proj > wind) 
    message (prog, "projection cannot exceed data window", ERROR);

  if (argc == 1) {
    ifstream* inputfile = new ifstream (*argv);
    if (inputfile -> good()) {
      cin = *inputfile;
    } else {
      message (prog, "unable to open input file", ERROR);
    }
  }
}


static void getdata (istream&        file,
		     vector<double>& data,
		     const int       skip,
		     const int       wind)
// ---------------------------------------------------------------------------
// Read in all the data, with appropriate skips; check length at least
// fills data window.
// ---------------------------------------------------------------------------
{
  Stack<double> datastk;
  double        datum;
  int           i, Nt;

  i = 0; while (file >> datum) if (!(i++ % skip)) datastk.push (datum);

  if ((Nt = datastk.depth()) < wind) { 
    cerr << prog << ": insufficient data to fill window" << endl;
    exit (EXIT_FAILURE);
  }

  data.setSize (Nt);
  for (i = 0; i < Nt; i++) data[Nt-i-1] = datastk.pop();
}


static void covary  (vector<double>& data, 
		     vector<double>& cov ,
		     const int       n   )
// ---------------------------------------------------------------------------
// Construct covariance matrix from lagged data vectors.  Covariance
// matrix is symmetric: we construct it in Lapack packed format.
// ---------------------------------------------------------------------------
{
  const int    N = data.getSize() - n + 1;
  const int    P = n+((n*(n-1))>>1);
  register int i, j, k;

  cov.setSize (P);

  cov = 0.0;
  for (i = 0; i < n; i++)
    for (j = i; j < n; j++)
      for (k = 0; k < N; k++)
	cov (Lapack::pack_addr (i, j)) += data (k + i) * data (k + j);

  Blas::scal (P, 1.0/N, cov(), 1);
}


static void eigensys (vector<double>& cov ,
		      vector<double>& esys,
		      const int       n   )
// ---------------------------------------------------------------------------
// Compute eigenvalues/vectors of covariance matrix, sort into
// descending eigenvalue order. In esys, the first n locations contain
// eigenvalues, followed by n lots of n-long (orthonomal)
// eigenvectors.
// ---------------------------------------------------------------------------
{
  const int      No2 = n >> 1;
  int            i, b;
  vector<double> work (3 * n);
  double         tmp;

  esys.setSize ((n+1)*n);

  // -- Compute eigenvalues (in ascending order), with matching eigenvectors.

  Lapack::spev ("V", "U", n, cov(), esys(), esys() + n, n, work(), i);

  if (i) message (prog, "eigensystem solution failed", ERROR);

  // -- Reverse the order supplied by spev.

  for (i = 0; i < No2; i++) {
    b = n - i - 1;
    tmp     = esys[i];
    esys[i] = esys[b];
    esys[b] = tmp;
    Veclib::copy (n, esys() + n * (i + 1), 1, work(),               1);
    Veclib::copy (n, esys() + n * (b + 1), 1, esys() + n * (i + 1), 1);
    Veclib::copy (n, work(),               1, esys() + n * (b + 1), 1);
  }
}


static void putsys (const vector<double>& esys,
		    const int             n    )
// ---------------------------------------------------------------------------
// Open a file called eigen.sys, print up eigensystem in it.
// 
// First line of file gives number of eigenpairs.  Then for each
// eigenpair we have a blank line, eigenvalue, followed by
// corresponding eigenvector.
// ---------------------------------------------------------------------------
{
  int       i, j;
  ofstream  file ("eigen.sys");

  file << n << endl;

  for (i = 0; i < n; i++) {
    file << endl << esys (i) << endl;
    for (j = 0; j < n; j++)
      file << esys (n + i*n + j) << endl;
  }
  
  file.close();
}


static void project (const vector<double>& data,
		     const int             p   ,
		     const vector<double>& esys, 
		     const int             n   )
// ---------------------------------------------------------------------------
// Project lagged data onto p leading eigenvectors of covariance
// matrix, print up. Each row gives the projection onto the first p
// vectors, and there are N rows.
// ---------------------------------------------------------------------------
{
  const int     N = data.getSize() - n + 1;
  const double* vect = esys() + n;
  int           i, j;

  if (!p) return;

  for (i = 0; i < N; i++) {
    for (j = 0; j < p; j++) {
      cout << "  " <<  Blas::dot (n, data() + i, 1, vect + n * j, 1);
    }
    cout << endl;
  }
}

