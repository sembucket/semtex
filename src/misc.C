///////////////////////////////////////////////////////////////////////////////
//misc.C: miscellaneous routines for I/O, memory management.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>


ostream& printVector (ostream&     strm,
		      const char*  fmt , 
		      const int    ntot,
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
  char    routine[] = "printVector";
  int     nvect;
  va_list ap;

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


char* upperCase (char *s)
// ---------------------------------------------------------------------------
// Uppercase characters in string.
// ---------------------------------------------------------------------------
{
  char *z(s); while (*z = toupper (*z)) z++; return s;
}
