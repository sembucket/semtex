/*****************************************************************************
 * Domain.C:  implement Domain class.
 *****************************************************************************/

// $Id$


#include "Fem.h"


static istream&   readOptions (istream&);
static istream&   readIparams (istream&);
static istream&   readFparams (istream&);





Domain::Domain (ifstream& file, const char* session)
// ---------------------------------------------------------------------------
// Construct a new Domain with a single Field from a named file.
//
// The named parts of input (separated by blank lines) are, in order:
//   ** OPTION               parameters
//   ** INTEGER              parameters
//   ** FLOATing point       parameters
//   ** BOUNDARY CONDITIONs
//   ** CURVED   EDGEs
//   ** MESH     INFORMATION
//
// On return, all information in Element and Boundary lists is set for u[0].
// ---------------------------------------------------------------------------
{
  char routine[] = "Domain::Domain";
  char s[StrMax];

  domain_name = new char[StrMax];
  strcpy (domain_name, session);

  domain_step = 0;
  domain_time = 0.0;
  nfield      = 1;

  u    = new Field* [nfield];
  u[0] = new Field;

  // -- Scan the parts of a session file.
   
  nextBlock (file, s);

  if (strstr (s, "**") && strstr (s, "OPTION"))    {
    readOptions (file);     
    nextBlock   (file, s); 
  }
  
  if (strstr (s, "**") && strstr (s, "INTEGER"))   {
    readIparams (file);     
    nextBlock   (file, s);
  }
  
  if (strstr (s, "**") && strstr (s, "FLOAT"))     {
    readFparams (file);     
    nextBlock   (file, s); 
  }
 
  if (strstr (s, "**") && strstr (s, "BOUNDARY") && strstr (s, "CONDITION")) { 
    if (!BCmanager::read (file)) message (routine, "no BCs set", ERROR);
    nextBlock       (file, s); 
  } else
    message (routine, "can't find boundary conditions", ERROR);
   
  if (strstr (s, "**") && strstr (s, "CURVED") && strstr (s, "EDGE"))    {
    // readCurvedEdges (file); 
    nextBlock       (file, s); 
  }
  
  if (strstr (s, "**") && strstr (s, "MESH") && strstr (s, "INFORMATION")) {
    ifstream* meshfile = altFile (file);
    u[0] -> readMesh (*meshfile, iparam ("N_POLY"));

    if (!u[0] -> nEl()) message (routine, "no element information set", ERROR);
    if (*meshfile != file) { meshfile -> close (); delete meshfile; }
  } else
    message (routine, "can't find mesh information", ERROR);

  file.close();

  // -- Complete computation of mesh geometric information.

  u[0] -> mapElements ();

  // -- Compute domain boundary edge information.
  
  u[0] -> buildBoundaries ();

  // -- Input mesh topology information, do RCM renumbering.

  u[0] -> readConnect (session);
  u[0] -> renumber    ();

  // -- Clean up.

  state_file.open   (strcat (strcpy (s, domain_name), ".sta"));
  history_file.open (strcat (strcpy (s, domain_name), ".his"));
  field_file.open   (strcat (strcpy (s, domain_name), ".fld"));

  if (!state_file)   message (routine, "can't open state file",   ERROR);
  if (!history_file) message (routine, "can't open history file", ERROR);
  if (!field_file)   message (routine, "can't open field file",   ERROR);
}





void Domain::addField (Field* F)
// ---------------------------------------------------------------------------
// Add a new Field pointer to Domain's array.
// ---------------------------------------------------------------------------
{
  Field** X = new Field* [nfield + 1];

  for (int i(0); i < nfield; i++) X[i] = u[i];
  X[nfield] = F;

  delete [] u;
  u = X;

  nfield++;
}





void Domain::restart ()
// ---------------------------------------------------------------------------
// If a restart file "name".rst can be found, use it to load
// velocity fields.  Otherwise initialize velocity fields to zero.
// ---------------------------------------------------------------------------
{
  char  restartfile[StrMax];
  strcat (strcpy (restartfile, domain_name), ".rst");

  ifstream file (restartfile);

  if   (file) file >> *this;
  else        for (int i(0); i < iparam ("N_VAR"); i++) *u[i] = 0.0;
}





static istream& readOptions (istream& istr)
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





static istream& readIparams (istream& istr)
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





static istream& readFparams (istream& istr)
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
