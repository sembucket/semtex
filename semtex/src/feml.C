///////////////////////////////////////////////////////////////////////////////
// feml.C:  Finite Element Markup Language (FEML) routines.
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
// 2. There are reserved keywords; see Feml.h (and below).
// 3. There is a maximum number of keywords (KEYWORDS_MAX) set in Feml.h.
// 4. No restriction is placed on options, or on the contents of a section.
// 5. After seeking, input stream is repositioned just following keyword.
//
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <Utility.h>
#include <Feml.h>
#include <Femlib.h>


FEML::FEML (const char* name)
// ---------------------------------------------------------------------------
// Attempt to open named file, prescan to locate sections.
// Look for name and name.pdf.  Set pdf_root.  Load tokens.
// ---------------------------------------------------------------------------
{
  char  c, routine[] = "FEML::open";
  char  err[StrMax], key[StrMax], yek[StrMax];
  char* u;
  int   i, OK, N, found;

  char* reserved[] = {
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
    0
  };
  
  feml_file.open (name);

  if (!feml_file) {
    strcat (strcpy (err, name), ".feml");
    feml_file.open (err);
  }

  if (!feml_file) {
    sprintf (err, "couldn't open file %s or %s.feml", name, name);
    message (routine, err, ERROR);
  }

  feml_root = strdup (name);
  if (strstr (name, ".feml")) feml_root[strlen (name) - 5] = '\0';

  for (i = 0; reserved[i] && i < KEYWORD_MAX; i++) {
    keyWord[i] = strdup (reserved[i]);
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

      for (found = 0, i = 0; !found && keyWord[i]; i++) {

	if (strcmp (key, keyWord[i]) == 0) {

	  // -- Locate closing '>'.

	  keyPosn[i] = feml_file.tellg ();
	  while (!found && feml_file >> c) found = c == '>';

	  if (!found) {
	    sprintf (err, "closing '>' not found for keyword %s", key);
	    message (routine, err, ERROR);
	  }

	  // -- Locate closing "</key>".

	  OK = 0;
      
	  while (!OK && feml_file >> c) {
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


int FEML::seek (const char* keyword)
// ---------------------------------------------------------------------------
// Look for keyword in stored table.
// If present, stream is positioned after keyword and 1 is returned.
// If not, stream is rewound and 0 is returned.
// ---------------------------------------------------------------------------
{
  register int i, found = 0;

  for (i = 0; !found && keyWord[i]; i++)
    found = (strstr (keyword, keyWord[i]) != 0 &&
	                      keyPosn[i]  != 0);

  if   (!found) {
    feml_file.clear ();  
    feml_file.seekg (0);
    return 0;
  } else{
    feml_file.clear ();
    feml_file.seekg (keyPosn[i - 1]);
  }

  return 1;
}


int FEML::attribute (const char* tag ,
		     const char* attr)
// ---------------------------------------------------------------------------
// Tag attributes are given as options in form <tag attr=int [attr=int ...]>
// Return integer value following '='.  No whitespace allowed in attributes.
// On return, FEML's stream is set to start of next line.
// ---------------------------------------------------------------------------
{
  char  routine[] = "FEML::attribute", buf[StrMax], err[StrMax];
  char* v;
  int   n = 0;

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

  feml_file.ignore (StrMax, '\n');

  return n;
}


int FEML::tokens ()
// ---------------------------------------------------------------------------
// Initialize femlib parser and install token table.
// Return 0 if no TOKEN section is found.
// NUMBER attribute ignored if present.
// ---------------------------------------------------------------------------
{
  char           buf[StrMax];
  register char* u;

  Femlib::prep();
 
  if (seek ("TOKENS")) {
    feml_file.ignore (StrMax, '\n');

    while (feml_file.getline (buf, StrMax)) {
      if (strstr (buf, "=")) Femlib::value (buf);
      u = buf; while (*u = toupper (*u)) u++;
      if (strstr (buf, "TOKENS")) break;
    }
    
    return 1;
  }
  
  return 0;
}
