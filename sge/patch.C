///////////////////////////////////////////////////////////////////////////////
// mixpatch.C: implement class functions for patching Mixed BCs on two domains.
//
// Copyright (C) 1999 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "Sem.h"


MixPatch::MixPatch (const Domain* d1,
		    const Domain* d2) :
// ---------------------------------------------------------------------------
// Build storage for patch BCs.
// ---------------------------------------------------------------------------
  np (Geometry::nP()),
  nz (Geometry::nZProc())
{
  const char routine[] = "MixPatch::MixPatch";
  char                    err[StrMax];
  const BoundarySys*      bsys;
  const Boundary*         B;
  Patch*                  P;
  vector<real>            work (np + np);
  real                    *tx = work(), *ty = tx + np;
  integer                 i, j, bad, nsurf, found;

  // -- Build storage by traversing the BCs of d1, adding for each mixed BC.

  bsys  = d1 -> b[0];
  nsurf = bsys -> nSurf();
  const vector<Boundary*>& BC1 = bsys -> BCs();

  for (i = 0; i < nsurf; i++) {
    B = BC1[i];
    if (strstr (B -> group(), "mixed")) {
      P         = new Patch;
      P -> b1   = B;
      P -> b2   = 0;
      P -> x    = new real [2 * np + 2 * np * nz];
      P -> y    = P -> x    + np;
      P -> val1 = P -> y    + np;
      P -> val2 = P -> val1 + np * nz;
      
      B -> geometry (P -> x, P -> y);
      Veclib::zero (2 * np * nz, P -> val1, 1);

      patches.add (P);
    }
  }

  ListIterator<Patch*> p (patches);

  // -- Traverse BCS of d2 and attempt to match up with those found
  // with d1.  The match is taken on spatial location of boundary
  // nodes.  We are assured that directions of traverse of patches are
  // reversed in d1 & d2's Boundaries.

  bsys  = d2 -> b[0];
  nsurf = bsys -> nSurf();
  const vector<Boundary*>& BC2 = bsys -> BCs();

  for (i = 0; i < nsurf; i++) {
    B = BC2[i];
    if (strstr (B -> group(), "mixed")) {
      found = 0;
      for (p.reset(); !found && p.more(); p.next()) {
	P = p.current();
	B -> geometry (tx, ty);
	for (bad = 0, j = 0; !bad && j < np; j++)
	  bad = ( fabs (P -> x[np - j - 1] - tx[j]) > EPSDP ||
		  fabs (P -> y[np - j - 1] - ty[j]) > EPSDP );
	if (!bad) { found = 1; break; }
      }

      if (found)
	P -> b2 = B;
      else {
	sprintf (err, "couldn't match %s's mixed BC %1d",
		 d2 -> name, B -> ID() + 1);
	message (routine, err, ERROR);
      }
    }
  }

  // -- Check that all patches have two Boundaries.

  for (p.reset(); p.more(); p.next())
    if (!(p.current() -> b1 && p.current() -> b2))
      message (routine, "mixed BC areas do not meet", ERROR);
}

