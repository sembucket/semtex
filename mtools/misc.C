///////////////////////////////////////////////////////////////////////////////
// misc.C:  Miscellaneous I/O routines.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <qmesh.h>


void error (const char* routine,
	    const char* text   ,
	    const lev&  level  )
// ---------------------------------------------------------------------------
// Error message handler, but without calls to graphics routines.
// ---------------------------------------------------------------------------
{
  switch (level) {
  case WARNING:
    cerr << "WARNING: " << routine << ": " << text << endl;
    return;
  case REMARK:
    cerr << text << endl;
    return;
  case ERROR:
    cerr << "ERROR: " << routine << ": " << text << endl;
    exit (EXIT_FAILURE);
    break;
  }
}


istream& seekBlock (istream&    strm,
		    const char* name)
// ---------------------------------------------------------------------------
// Search input file stream for a block starting with name (case insensitive).
// Then read on until first "{" is encountered.
// ---------------------------------------------------------------------------
{
  char  routine[] = "seekBlock";
  char  s[StrMax], uname[StrMax];
  int   c;

  strcpy    (uname, name);
  upperCase (uname);

  strm.seekg(0).clear();

  while (strm >> s) {
    upperCase (s);
    if (strcmp (s, uname) == 0) {
      while ((c = strm.get ()) != '{' && c != '}' && c != EOF);
      if (c == EOF) 
	error (routine, "reached EOF while looking for '{'", ERROR);
      if (c == '}')
	error (routine,   "found '}' while looking for '{'", ERROR);
      return strm;
    } else
      strm.getline (s, StrMax);
  }

  sprintf (s, "failed to locate block \"%s\"", name);
  error   (routine, s, ERROR);
  return  strm;
}


istream& endBlock (istream& strm)
// ---------------------------------------------------------------------------
// Advance to next "}" in file.
// ---------------------------------------------------------------------------
{
  char  routine[] = "endBlock";
  int   c;

  while ((c = strm.get ()) != '}' && c != EOF);
  if (c == EOF) error (routine, "reached EOF while looking for '}'", ERROR);
  
  return strm;
}


char* upperCase (char *s)
// ---------------------------------------------------------------------------
// Uppercase characters in string.
// ---------------------------------------------------------------------------
{
  char *z(s); while (*z = toupper (*z)) z++; return s;
}


void message (const char* routine,
	      const char* text   ,
	      const lev&  level  )
// ---------------------------------------------------------------------------
// Error message handler for all modules that may run graphics commands.
// ---------------------------------------------------------------------------
{
  switch (level) {
  case WARNING:
    cerr << "WARNING: " << routine << ": " << text << endl;
    return;
  case REMARK:
    cerr << text << endl;
    return;
  case ERROR:
    cerr << "ERROR: " << routine << ": " << text << endl;
#if defined(GRAPHICS)
    if (graphics) stopGraphics ();
#endif
    exit (EXIT_FAILURE);
    break;
  }
}
