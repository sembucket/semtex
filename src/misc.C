/*****************************************************************************
 * misc.C: miscellaneous routines for I/O, memory management.
 *****************************************************************************/

static char
RCSid[] = "$Id$";

#include "Fem.h"
#include <time.h>


ostream& printVector (ostream&     strm,
		      const char*  fmt , 
		      int          ntot,
		                   ... )
// ---------------------------------------------------------------------------
// Print up a variable number of numeric vectors on strm, in columns.
//
// The format specifier is gives the number and type of the vectors, with
// type specified by the first character.  The vectors must all be of the
// same type & length.
//
// Allowed types are: int ("i"), real ("r").
// Examples: four integer vectors ==> fmt is "iiii".  Two reals ==> "rr".
// 
// Vectors are printed in a fixed field width of 15, regardless of type.
// ---------------------------------------------------------------------------
{
  char     routine[] = "printVector";
  int      nvect;
  va_list  ap;

  nvect = strlen (fmt);
  if (! nvect   ) message (routine, "empty format string",   ERROR  );
  if (nvect > 10) message (routine, "more than 10 vectors?", WARNING);
  if (ntot  <  0) message (routine, "ntot < 0",              ERROR  );

  switch (fmt[0]) {

  case 'i': {
    int** u = new int* [nvect];
    va_start (ap, ntot);
    for (int k = 0; k < nvect; k++) u[k] = va_arg (ap, int*);
    va_end (ap);
    for (register int l = 0; l < ntot; l++) {
      for (register int j = 0; j < nvect; j++)
	strm << setw(15) << u[j][l];
      strm << endl;
    }
    delete [] u;
    break;
  }
  case 'r': {
    real** u = new real* [nvect];
    va_start (ap, ntot);
    for (int k = 0; k < nvect; k++) u[k] = va_arg (ap, real*);
    va_end (ap);
    for (register int l = 0; l < ntot; l++) {
      for (register int j = 0; j < nvect; j++)
	strm << setw(15) << u[j][l];
      strm << endl;
    }
    delete [] u;
    break;
  }
  default:
    message (routine, fmt, ERROR);
    break;
  }

  if (!strm) message (routine, "output failed", ERROR);

  return strm;
}


ifstream&  seekBlock (ifstream& strm, const char* name)
// ---------------------------------------------------------------------------
// Search input file stream for a block starting with name (case insensitive).
// Then read on until first "{" is encountered.
// ---------------------------------------------------------------------------
{
  char  routine[] = "seekBlock";
  char  s[StrMax], uname[StrMax];
  int   c;

  strcpy     (uname, name);
  upperCase  (uname);
  strm.seekg (0, ios::beg);

  while (strm >> s) {
    upperCase (s);
    if (strcmp (s, uname) == 0) {
      while ((c = strm.get ()) != '{' && c != '}' && c != EOF);
      if (c == EOF) 
	message (routine, "reached EOF while looking for '{'", ERROR);
      if (c == '}')
	message (routine,   "found '}' while looking for '{'", ERROR);
      return strm;
    } else
      strm.getline (s, StrMax);
  }

  sprintf (s, "failed to locate block \"%s\"", name);
  message (routine, s, ERROR);
  return  strm;
}


istream&  endBlock (istream& strm)
// ---------------------------------------------------------------------------
// Advance to next "}" in file.
// ---------------------------------------------------------------------------
{
  char  routine[] = "endBlock";
  int   c;

  while ((c = strm.get ()) != '}' && c != EOF);
  if (c == EOF) message (routine, "reached EOF while looking for '}'", ERROR);
  
  return strm;
}


char*  upperCase (char *s)
// ---------------------------------------------------------------------------
// Uppercase characters in string.
// ---------------------------------------------------------------------------
{
  char *z(s); while (*z = toupper (*z)) z++; return s;
}


ifstream*  altFile (ifstream& ist)
// ---------------------------------------------------------------------------
// Look in file for a new filename.  Open & return a pointer if it exists.
// ---------------------------------------------------------------------------
{
  char routine[] = "altFile";
  char s1[StrMax];

  ist.getline(s1, StrMax);

  if (   strstr (s1, "FILE") 
      || strstr (s1, "File") 
      || strstr (s1, "file") ) {

    ifstream  *newfile;
    char       s2[StrMax];

    if (!(sscanf (s1, "%*s, %s", s2))) {
      sprintf (s2, "couldn't get element file name from string: %s", s1);
      message (routine, s2, ERROR);
    } else if (!(newfile = new ifstream (s2))) {
      sprintf (s1, "couldn't find alternate file %s", s2);
      message (routine, s1, ERROR);
    }

    return newfile;

  } else
    return &ist;
}


int*  ivector (long len)
// ----------------------------------------------------------------------------
// Return zero-offset vector of int.
// ----------------------------------------------------------------------------
{
  int* v = new int [len];
  
  return v;
}


real*  rvector (long len)
// ----------------------------------------------------------------------------
// Return zero-offset vector of real.
// ----------------------------------------------------------------------------
{
  real* v = new real [len];
  
  return v;
}


void freeVector (int* v)
// ----------------------------------------------------------------------------
// Free a zero-offset vector.
// ----------------------------------------------------------------------------
{
  delete [] v;
}


void freeVector (real* v)
// ----------------------------------------------------------------------------
// Free a zero-offset vector.
// ----------------------------------------------------------------------------
{
  delete [] v;
}


int** imatrix (long nrow, long ncol)
// ---------------------------------------------------------------------------
// Allocate a 2D row-major, zero-offset, contiguous matrix of int.
// ---------------------------------------------------------------------------
{
  int** m = new int* [nrow];
  m[0]    = new int  [nrow*ncol];
  for (register long i = 1; i < nrow; i++) m[i] = m[i - 1] + ncol;

  return m;
}


real** rmatrix (long nrow, long ncol)
// ---------------------------------------------------------------------------
// Allocate a 2D row-major, zero-offset, contiguous matrix of real.
// ---------------------------------------------------------------------------
{
  real** m = new real* [nrow];
  m[0]     = new real  [nrow*ncol];
  for (register long i = 1; i < nrow; i++) m[i] = m[i - 1] + ncol;

  return m;
}


void freeMatrix (int** m)
// ---------------------------------------------------------------------------
// Release storage for int matrix.
// ---------------------------------------------------------------------------
{
  delete [] m[0];
  delete [] m;
}


void freeMatrix (real** m)
// ---------------------------------------------------------------------------
// Release storage for real matrix.
// ---------------------------------------------------------------------------
{
  delete [] m[0];
  delete [] m;
}


istream& readOptions (istream& istr)
// ---------------------------------------------------------------------------
// Read solution options from istr, set into external list.
//
// Option parameters may be either integer or character (the first char of
// a string), but in either case the integer representation is stored.  This
// means that string option parameters can be distinguished if the first
// character is unique (case sensitive).
// ---------------------------------------------------------------------------
{
  char   routine[] = "readOptions";
  char   p[StrMax], s[StrMax], err[StrMax];
  int    n, i;
  
  istr >> s;
  istrstream (s, strlen (s)) >> n;
  if (n < 0) {
    sprintf (err, "expected integer number of options, read: %s", s);
    message (routine, err, ERROR);
  }
  if (n > 20) {
    sprintf (err, "found a large number of options (%1d) in: %s", s);
    message (routine, err, WARNING);
  }
  istr.getline (s, StrMax);
  upperCase    (s);
  if (!strstr (s, "OPTION")) {
    sprintf (err, "can't locate 'OPTION' string in: %s", s);
    message (routine, err, ERROR);
  }

  while (n--) {
    istr >> p;
    istr >> s;
    if   (isalpha (p[0])) i = (int) p[0];
    else                  i = atoi (p);
    setOption (s, i);
  }

  return istr;
}


istream& readIparams (istream& istr)
// ---------------------------------------------------------------------------
// Read integer parameters from istr, set into external list.
// ---------------------------------------------------------------------------
{
  char   routine[] = "readIparams";
  char   s[StrMax], err[StrMax];
  int    n, i;
  
  istr >> s;
  istrstream (s, strlen (s)) >> n;
  if (n < 0) {
    sprintf (err, "expected number of parameters, read: %s", s);
    message (routine, err, ERROR);
  }
  if (n > 20) {
    sprintf (err, "found a large number of parameters (%1d) in: %s", s);
    message (routine, err, WARNING);
  }
  istr.getline (s, StrMax);
  upperCase    (s);
  if (!strstr (s, "INTEGER")) {
    sprintf (err, "can't locate 'INTEGER' string in: %s", s);
    message (routine, err, ERROR);
  }

  while (n--) {
    istr >> i;
    istr >> s;
    setIparam (s, i);
  }

  return istr;
}


istream& readFparams (istream& istr)
// ---------------------------------------------------------------------------
// Read floating-point parameters from istr, set into external list.
//
// Values may be defined in terms of previously-defined symbols.
// Symbols are case-sensitive.
//
// Example:
// 500      Re
// 1/Re     KINVIS  (sets dparam "KINVIS" = 0.002).
// ---------------------------------------------------------------------------
{
  char   routine[] = "readFparams";
  char   p[StrMax], s[StrMax], err[StrMax];
  int    n;
  
  istr >> s;
  istrstream (s, strlen (s)) >> n;
  if (n < 0) {
    sprintf (err, "expected number of parameters, read: %s", s);
    message (routine, err, ERROR);
  }
  if (n > 20) {
    sprintf (err, "found a large number of parameters (%1d) in: %s", s);
    message (routine, err, WARNING);
  }
  istr.getline (s, StrMax);
  upperCase    (s);
  if (!strstr (s, "FLOAT")) {
    sprintf (err, "can't locate 'FLOAT' string in: %s", s);
    message (routine, err, ERROR);
  }

  while (n--) {
    istr >> p;
    istr >> s;
    setDparam (s, interpret (p));
  }

  return istr;
}


Mesh*  preProcess (ifstream& strm)
// ---------------------------------------------------------------------------
// Get run-time parameters, BCs & Mesh from strm.
//
// Set default options, parameters, according to problem.
// ---------------------------------------------------------------------------
{
  char   routine[] = "preProcess";
  int    verbose   = option ("VERBOSE");
  Mesh*  M         = new Mesh;

  seekBlock   (strm, "parameter");
  readOptions (strm);     
  readIparams (strm);     
  readFparams (strm);     
  endBlock    (strm);
  
  if (!iparam ("N_POLY")) {
    setIparam ("N_POLY", 5);
    message (routine, "setting default N_POLY = 5", WARNING);
  }

  if (option ("PROBLEM") == NAVIERSTOKES) {
    if (!iparam ("N_STEP")) {
      setIparam ("N_STEP", 10);
      message (routine, "setting default N_STEP = 10", WARNING);
    }   
    if (!iparam ("IO_FLD")) {
      setIparam ("IO_FLD", iparam ("N_STEP"));
      message (routine, "setting default IO_FLD = N_STEP", WARNING);
    }
    if (!iparam ("N_TIME")) {
      setIparam ("N_TIME", 1);
      message (routine, "setting default N_TIME = 1", WARNING);
    }
    if (dparam ("DELTAT") == 0.0) {
      setDparam ("DELTAT", 0.01);
      message (routine, "setting default DELTAT = 0.01", WARNING);
    }
    if (dparam ("KINVIS") == 0.0) {
      setDparam ("KINVIS", 1.0);
      setDparam ("DELTAT", 0.01);
      message (routine, "setting default KINVIS = 1.0", WARNING);
    }
  } else if (option ("PROBLEM") == HELMHOLTZ) {
    setIparam ("N_STEP", 1);
    setIparam ("IO_FLD", 1);
    setDparam ("DELTAT", 0.0);
  }

  if (verbose) {
    message (routine, " -- OPTION PARAMETERS:",         REMARK); showOption ();
    message (routine, " -- INTEGER PARAMETERS:",        REMARK); showIparam ();
    message (routine, " -- FLOATING POINT PARAMETERS:", REMARK); showDparam ();
  }

  seekBlock (strm, "boundary");
  if (!BCmanager::read (strm)) message (routine, "no BCs set", ERROR);
  endBlock  (strm); 

  if (verbose) {
    message (routine, " -- BOUNDARY CONDITIONS:", REMARK);
    BCmanager::print ();
  }

  seekBlock (strm, "mesh");
  strm >> *M;
  endBlock  (strm);

  if (verbose > 1) {
    message (routine, " -- MESH ASSEMBLY INFORMATION:", REMARK);
    Mesh::printAssembly (*M);
  }

  return M;
}


