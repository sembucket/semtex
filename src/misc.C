/*****************************************************************************
 * misc.C: miscellaneous routines for I/O, memory management.
 *****************************************************************************/

// $Id$

#include "Fem.h"
#include <time.h>





ostream& printVector (ostream&     strm,
		      const char*  fmt , 
		      int          ntot,
		                   ... )
// ---------------------------------------------------------------------------
// Print up a variable number of vectors on strm, in columns.
//
// The format specifier is gives the number and type of the vectors, with
// type specified by the first character.
//
// Allowed types are: int ("i"), real ("r").
// Examples: four integer vectors ==> fmt is "iiii".  Two reals ==> "rr".
// 
// Vectors are printed in a fixed field width of 15, regardless of type.
// ---------------------------------------------------------------------------
{
  char      routine[] = "printVector";
  int       nvect;
  va_list   ap;

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





istream&  nextBlock (istream& strm, char* s)
// ---------------------------------------------------------------------------
// Advance to start of the next block of information in file, skipping
// empty lines and comment lines (lines starting with '#').
// Then uppercase the new start line.
// ---------------------------------------------------------------------------
{
  while (strm.getline (s, StrMax))
    if (s[0] != '#') break;
  upperCase (s);

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





int  quadComplete (int dim, int np)
// ---------------------------------------------------------------------------
// Return the number of Gauss-Legendre quadrature points sufficient to
// achieve the full rate of convergence for tensor-product element bases.
//
// Dim is the number of space dimensions, np the number of points defining
// basis polynomials.
//
// References: Hughes \S 4.1, Strang & Fix \S 4.3.
// ---------------------------------------------------------------------------
{
  int  n, ktot;
  

  ktot = (dim + 1)*(np - 1) - 2;
  n = (ktot & 0) ? ktot + 2 : ktot + 1;
  n >>= 1;

  return max (n, 2);
}





ostream& operator << (ostream& strm, Domain& D)
// ---------------------------------------------------------------------------
// Output all Domain field variables on ostream in prism-compatible form.
// A friend of Field.
// ---------------------------------------------------------------------------
{
  char *hdr_fmt[] = { 
    "%-25s "    "Session\n",
    "%-25s "    "Created\n",
    "%-25s "    "Nr, Ns, Nz, Elements\n",
    "%-25d "    "Step\n",
    "%-25.6g "  "Time\n",
    "%-25.6g "  "Time step\n",
    "%-25.6g "  "Kinvis\n",
    "%-25.6g "  "Beta\n",
    "%-25s "    "Fields Written\n",
    "%-25s "    "Format\n"
  };

  char      routine[] = "Domain << operator";
  char      s1[StrMax], s2[StrMax];
  int       np   (D.u[0] -> element_list.first() -> nKnot());
  int       ntot (D.u[0] -> nTot());
  time_t    tp   (::time(0));

  sprintf (s1, hdr_fmt[0], D.name());
  strm << s1;

  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  sprintf  (s1, hdr_fmt[1], s2);
  strm << s1;

  sprintf (s2, "%1d %1d %1d %1d", np, np, 1, D.u[0] -> nEl());
  sprintf (s1, hdr_fmt[2], s2);
  strm << s1;

  sprintf (s1, hdr_fmt[3], D.step());
  strm << s1;

  sprintf (s1, hdr_fmt[4], D.time());
  strm << s1;

  sprintf (s1, hdr_fmt[5], dparam ("DELTAT"));
  strm << s1;

  sprintf (s1, hdr_fmt[6], dparam ("KINVIS"));
  strm << s1;

  sprintf (s1, hdr_fmt[7], dparam ("BETA"));
  strm << s1;

  int k;
  for (k = 0; k < D.nField(); k++) s2[k] = D.u[k] -> name();
  s2[k] = '\0';
  sprintf (s1, hdr_fmt[8], s2);
  strm << s1;

  sprintf (s1, hdr_fmt[9], (option ("BINARY")) ? "binary" : "ASCII");
  strm << s1;

  if (option ("BINARY")) {
    for (register int n(0); n < D.nField(); n++)
      strm.write((char*) D.u[n]->data, ntot * sizeof (real));
  } else {
    strm.setf (ios::scientific, ios::floatfield); strm.setf (ios::uppercase);
    for (register int i(0); i < ntot; i++) {
      for (register int n(0); n < D.nField(); n++)
        strm << setw (14) << D.u[n] -> data[i];
      strm << endl;
    }
  }

  if (!strm) message (routine, "failed writing field file", ERROR);
  strm << flush;

  return strm;
}





istream& operator >> (istream& strm, Domain& D)
// ---------------------------------------------------------------------------
// Input all Domain field variables from prism-compatible istream.
// A friend of Field.
// ---------------------------------------------------------------------------
{
  char          routine[] = "strm >> Domain";
  int           i, np, ns, nz, nel, ntot;
  char          s[StrMax];

  for (i = 0; i < 3; i++) strm.getline (s, StrMax);

  sscanf (s, "%d %d %d %d", &np, &ns, &nz, &nel);
  if ((np != D.u[0] -> element_list.first() -> nKnot()) ||
      (ns != D.u[0] -> element_list.first() -> nKnot()))
    message (routine, "element size mismatch",       ERROR);
  if (nz != 1)
    message (routine, "number of z planes mismatch", ERROR);
  if (nel != D.u[0] -> nEl())
    message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != D.u[0] -> nTot ())
    message (routine, "declared sizes mismatch",     ERROR);
  
  for (i = 3; i < 10; i++) strm.getline (s, StrMax);
  
  if (strstr (s, "binary")) {
    for (i = 0; i < D.nField () - 1; i++) {
      strm.read ((char *) D.u[i] -> data, ntot * sizeof (real));
      if (!strm)
	message (routine, "unable to read field from strm", ERROR);
    }
    
  } else if (strstr (s, "ASCII")) {
    for (register int j = 0; j < ntot; j++) {
      strm.getline (s, StrMax);
      if (!strm)
	message (routine, "premature EOF", ERROR);
      for (register int n = 0; n < D.nField () - 1; n++)
	if (sscanf (s, "%lf", D.u[n] -> data + j) < 1)
	  message (routine, "unable to read field from strm", ERROR);      
    }
  }
    
  return strm;
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
  char   s1[StrMax], s2[StrMax];
  int    i;
  
  istr.getline (s1, StrMax).getline (s1, StrMax);

  do {
    if (isalpha (s1[0])) {
      i = (int) s1[0];
      if (sscanf (s1, "%*s %s", s2) != 1) {
	sprintf (s2, "unable to scan character parameter from string %s", s1);
	message (routine, s2, ERROR);
      }
    } else {
      if (sscanf (s1, "%d %s", &i, s2) != 2) {
	sprintf (s2, "unable to scan integer parameter from string %s", s1);
	message (routine, s2, ERROR);
      }
    }

    setOption (s2, i);
    
    istr.getline(s1, StrMax);
  } while (s1[0]);

  return istr;
}





istream& readIparams (istream& istr)
// ---------------------------------------------------------------------------
// Read integer parameters from istr, set into external list.
// ---------------------------------------------------------------------------
{
  char   routine[] = "readIparams";
  char   s1[StrMax], s2[StrMax];
  int    i;

  istr.getline(s1, StrMax).getline(s1, StrMax);

  do {
    if (sscanf (s1, "%d %s", &i, s2) != 2) {
      sprintf (s2, "unable to scan integer parameter from string %s", s1);
      message (routine, s2, ERROR);
    }
    setIparam (s2, i);
    
    istr.getline(s1, StrMax);
  } while (s1[0]);

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
  char   s1[StrMax], s2[StrMax], s3[StrMax];

  istr.getline(s1, StrMax).getline(s1, StrMax);

  do {
    if (sscanf (s1, "%s %s", s2, s3) != 2)
      message (routine, strcat(s1, "?"), ERROR);

    setDparam(s3, interpret (s2));
    
    istr.getline(s1, StrMax);
  } while (s1[0]);

  return istr;
}





Mesh*  preProcess (istream& strm)
// ---------------------------------------------------------------------------
// Get run-time parameters, BCs & Mesh from strm.
//
// The named parts of input (separated by blank lines) are, in order:
//   ** OPTION               parameters
//   ** INTEGER              parameters
//   ** FLOATing point       parameters
//   ** BOUNDARY CONDITIONs
//   ** MESH     INFORMATION
// ---------------------------------------------------------------------------
{
  char   routine[] = "preProcess";
  char   s[StrMax];
  Mesh*  M = new Mesh;

  // -- Scan the parts of a session file.
   
  nextBlock (strm, s);

  if (strstr (s, "**") && strstr (s, "OPTION"))    {
    readOptions (strm);     
    nextBlock   (strm, s); 
  }
  
  if (strstr (s, "**") && strstr (s, "INTEGER"))   {
    readIparams (strm);     
    nextBlock   (strm, s);
  }
  
  if (strstr (s, "**") && strstr (s, "FLOAT"))     {
    readFparams (strm);     
    nextBlock   (strm, s); 
  }
 
  if (strstr (s, "**") && strstr (s, "BOUNDARY") && strstr (s, "CONDITION")) { 
    if (!BCmanager::read (strm)) message (routine, "no BCs set", ERROR);
    nextBlock       (strm, s); 
  } else
    message (routine, "can't find boundary conditions", ERROR);

  if (strstr (s, "**") && strstr (s, "MESH") && strstr (s, "INFORMATION")) {
    strm >> *M;
  } else
    message (routine, "can't find mesh information", ERROR);

  return M;
}
