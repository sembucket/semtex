///////////////////////////////////////////////////////////////////////////////
// feml.C:  Finite Element Markup Language (FEML) routines.
//
// Copyright (c) 1994,2003 Hugh Blackburn
//
// After initialization, FEML files are prescanned to find locations of
// keywords.  These locations are stored, and reset by the seek function.
//
// Keywords are set in SGML section-block notation, with open and close tags.
//
// <keyword [options]>
//   .
//   .
//   .
// </keyword>
//
// 1. Keywords case-insensitive on input, converted to uppercase internally.
// 2. There are reserved keywords; see feml (and below).
// 3. There is a maximum number of keywords (KEYWORDS_MAX) set in feml.
// 4. No restriction is placed on options, or on the contents of a section.
// 5. After seeking, input stream is repositioned just following keyword.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>

#include <cfemdef>
#include <utility_h>
#include <feml_h>
#include <femlib_h>


FEML::FEML (const char* session)
// ---------------------------------------------------------------------------
// Attempt to open session file, prescan to locate sections.  Set
// pdf_root.  Load tokens.
// ---------------------------------------------------------------------------
{
  const char routine[] = "FEML::FEML";
  char       c, err[STR_MAX], key[STR_MAX], yek[STR_MAX];
  char*      u;
  integer    i, N;
  bool       OK, found;

  const char* reserved[] = {
    "TOKENS",
    "FIELDS",
    "GROUPS",
    "BCS",
    "NODES",
    "ELEMENTS",
    "SURFACES",
    "CURVES",
    "USER",
    "HISTORY",
    "CUT",
    0
  };
  
  feml_file.open (session);

  if (!feml_file) {
    sprintf (err, "couldn't open session file %s", session);
    message (routine, err, ERROR);
  }

  strcpy ((feml_root = new char [strlen (session) + 1]), session);

  for (i = 0; reserved[i] && i < KEYWORD_MAX; i++) {
    keyWord[i] = strcpy ((new char [strlen (reserved[i]) + 1]), reserved[i]);
    keyPosn[i] = 0;
  }
  keyWord[i] = 0;

  if (i == KEYWORD_MAX) {
    sprintf (err, "Number of reserved keywords exceeds table size (%1d)", i);
    message (routine, err, ERROR);
  }

  while (feml_file >> c) {

    if (c == '<') {

      // -- Next word may be a keyword. 
      //    First massage it to get correct form.
      
      feml_file >> key;


      N = strlen (key);
      if (key[N - 1] == '>') {
	c = '>';
	feml_file.putback (c);
	key[N - 1] = '\0';
	N--;
      }

      u = key; while (*u = toupper (*u)) u++;

      // -- Check if key is a keyword and if so:
      //    1. install file position immediately following keyword in table;
      //    2. move on to find closing tag.

      for (found = false, i = 0; !found && keyWord[i]; i++) {

	if (strcmp (key, keyWord[i]) == 0) {

	  // -- Locate closing '>'.

	  keyPosn[i] = feml_file.tellg ();
	  while (!found && feml_file >> c) found = c == '>';

	  if (!found) {
	    sprintf (err, "closing '>' not found for keyword %s", key);
	    message (routine, err, ERROR);
	  }

	  // -- Locate closing "</key>".

	  OK = false;
      
	  while ((!OK) && (feml_file >> c)) {
	    if (c == '<') {
	      feml_file >> c;
	  
	      if (c == '/') {
		feml_file >> yek;
		
		N = strlen (yek);
		if (yek[N - 1] == '>') {
		  c = '>';
		  feml_file.putback (c);
		  yek[N - 1] = '\0';
		  N--;
		}

		u = yek; while (*u = toupper (*u)) u++;
	    
		if (OK = strcmp (key, yek) == 0) {
		  while (feml_file >> c) if (c == '>') break;

		  if (c != '>') {
		    sprintf (err, "closing '>' not found for /%s", key);
		    message (routine, err, ERROR);
		  }
		}
	      }
	    }
	  }

	  if (!OK) {
	    sprintf (err, "couldn't locate </%s> to match <%s>", key, key);
	    message (routine, err, ERROR);
	  }
	}
      }
    }
  }
  
  if (!found) message (routine, "no keywords located", ERROR);

  feml_file.clear ();		// -- Reset EOF error condition.
  feml_file.seekg (0);		// -- And rewind.

  tokens ();			// -- Initialize Femlib parser.
}


integer FEML::seek (const char* keyword)
// ---------------------------------------------------------------------------
// Look for keyword in stored table.
// If present, stream is positioned after keyword and 1 is returned.
// If not, stream is rewound and 0 is returned.
// ---------------------------------------------------------------------------
{
  register integer i;
  bool     found = false;

  for (i = 0; !found && keyWord[i]; i++)
    found = (strstr (keyword, keyWord[i]) != 0 &&
	     static_cast<integer>(keyPosn[i]) != 0);

  if (!found) {
    feml_file.clear ();  
    feml_file.seekg (0);
    return 0;
  } else {
    feml_file.clear ();
    feml_file.seekg (keyPosn[i - 1]);
  }

  return 1;
}


integer FEML::attribute (const char* tag ,
			 const char* attr)
// ---------------------------------------------------------------------------
// Tag attributes are given as options in form <tag attr=int [attr=int ...]>
// Return integer value following '='.  No whitespace allowed in attributes.
// On return, FEML's stream is set to start of next line.
// ---------------------------------------------------------------------------
{
  const char routine[] = "FEML::attribute";
  char       buf[STR_MAX], err[STR_MAX];
  char*      v;
  integer    n = 0;

  if (!seek (tag)) {
    sprintf (err, "couldn't locate tag %s in feml file", tag);
    message (routine, err, ERROR);
  }

  while (feml_file >> buf)
    if (strstr (buf, attr)) {
      if (buf[strlen (buf) - 1] == '>') buf[strlen (buf) - 1] = '\0';
      v = buf;
      while (*v && *v++ != '=');
      if (!(*v)) {
	sprintf (err, "attribute syntax error in %s", buf);
	message (routine, err, ERROR);
      }

      n = atoi (v);
      break;
    } else if (strchr (buf, '>')) {
      sprintf (err, "%s not found in tag %s", attr, tag);
      message (routine, err, ERROR);
    }

  feml_file.ignore (STR_MAX, '\n');

  return n;
}


bool FEML::tokens ()
// ---------------------------------------------------------------------------
// Install token table.  Return false if no TOKEN section is found.
// NUMBER attribute ignored if present.  Fix any inconsistent values.
// Parser must have been initialized before entry.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "FEML::tokens";
  char           buf[STR_MAX];
  register char* u;

  if (seek ("TOKENS")) {
    feml_file.ignore (STR_MAX, '\n');

    while (feml_file.getline (buf, STR_MAX)) {
      if (strstr (buf, "=")) Femlib::value (buf);
      u = buf; while (*u = toupper (*u)) u++;
      if (strstr (buf, "TOKENS")) break;
    }
    
    if (Femlib::value ("IO_FLD") > Femlib::ivalue ("N_STEP"))
      Femlib::ivalue ("IO_FLD", Femlib::ivalue ("N_STEP"));

    if (Femlib::ivalue ("N_TIME") > 3) {
      message (routine, "N_TIME too large, reset to 3", WARNING);
      Femlib::ivalue ("N_TIME", 3);
    }
    
    return true;
  }
  
  return false;
}
