/*****************************************************************************
 * Domain.C:  implement Domain class.
 *****************************************************************************/

// $Id$


#include "Fem.h"


Domain::Domain (Mesh& M, const char* session, int np)
// ---------------------------------------------------------------------------
// Construct a new Domain with a single Field from a named file.
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
  u[0] = new Field (M, np);

  // -- Complete computation of mesh geometric information.

  u[0] -> mapElements ();

  // -- Compute domain boundary edge information.
  
  u[0] -> buildBoundaries (M);

  // -- Input mesh topology information, do RCM renumbering.

  u[0] -> connect  (M, np);
  u[0] -> renumber ();

#ifdef DEBUG
  Field::printConnect (u[0]);
#endif

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
