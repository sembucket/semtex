/*****************************************************************************
 * domain.C:  implement Domain class.
 *****************************************************************************/

static char RCSid[] = "$Id$";

#include "Fem.h"


Domain::Domain (Mesh& M, const char* session, int np)
// ---------------------------------------------------------------------------
// Construct a new Domain with a single Field, using geometry from M.
//
// On return, all information in Element and Boundary lists is set for u[0].
// ---------------------------------------------------------------------------
{
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

  // -- Generate mesh topology information, do RCM renumbering.

  u[0] -> connect (M, np);

}


void Domain::openFiles ()
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  char  routine[] = "Domain::openFiles";
  char  s[StrMax];

  state_file.open   (strcat (strcpy (s, domain_name), ".sta"));
  history_file.open (strcat (strcpy (s, domain_name), ".his"));
  if (option ("CHKPOINT"))
    field_file.open (strcat (strcpy (s, domain_name), ".chk"));
  else
    field_file.open (strcat (strcpy (s, domain_name), ".fld"));

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


void  Domain::dump ()
// ---------------------------------------------------------------------------
// Check if a field-file write is required, carry out.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Domain::dump";
  int   verbose   = option ("VERBOSE");

  if (domain_step % iparam ("IO_FLD")) return;

  if (option ("CHKPOINT")) {
    char s[StrMax], b[StrMax], c[StrMax];

    field_file.close ();
    strcat (strcpy (s, domain_name), ".chk");
    strcat (strcpy (b, domain_name), ".chk.bak");

    sprintf (c, "mv ./%s ./%s", s, b);
    system  (c);

    field_file.open (s);
    if (!field_file)
      message (routine, "failed to open checkpoint file", ERROR);
  }

  if (verbose) message (routine, ": writing field dump", REMARK);
  field_file << *this;
}


void Domain::cleanup ()
// ---------------------------------------------------------------------------
// Do final output of domain variables if required, rename files if needed.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Domain::dump";
  int   verbose   = option ("VERBOSE");

  if (iparam ("N_STEP") % iparam ("IO_FLD")) {
    if (verbose) message (routine, ": writing field dump", REMARK);
    field_file << *this;
  }
  field_file.close ();

  if (option ("CHKPOINT")) {
    char s[StrMax], b[StrMax], c[StrMax];

    strcat  (strcpy (s, domain_name), ".chk");
    strcat  (strcpy (b, domain_name), ".fld");
    sprintf (c, "mv ./%s ./%s", s, b);
    system  (c);
  }
}

