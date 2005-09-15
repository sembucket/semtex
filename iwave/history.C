///////////////////////////////////////////////////////////////////////////////
// history.C
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// Routines to provide history point information at x, y, z locations.
//
// History points are nominated in the session file, e.g.
// <HISTORY NUMBER=2>
// #    tag   x     y    z
//      1     0.0   0.0  0.0
//      2     0.3  -0.1 -2
// </HISTORY>
//
// The z locational information is ignored for 2D.  For 3D the periodic
// assumption is used in the z-direction in case the z location does not
// fall within the first spatial period.
//
// Output to a file called session.his.
//
// *** This version is modified for iwave, to allow history points to
// remain fixed in azimuth relative to frame that holds spinning
// container. This is signalled by TOKEN LDV_FIXED (also need
// OMEGA_1).
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <sem.h>

const Element* HistoryPoint::locate (const real_t        x   ,
				     const real_t        y   ,
				     vector<Element*>& Esys,
				     real_t&             r   ,
				     real_t&             s   )
// ---------------------------------------------------------------------------
// Static class member function which tries to locate x, y, point within
// an element E, and return its location in r, s coordinates within E.
// ---------------------------------------------------------------------------
{
  register int_t i;
  const int_t    NEL   = Esys.size();
  const bool     guess = true;
  const Element* E;
  vector<real_t> work (max (2*Geometry::nTotElmt(), 5*Geometry::nP() + 6));

  for (E = 0, i = 0; i < NEL; i++)
    if (Esys[i] -> locate (x, y, r=0.0, s=0.0, &work[0], guess)) E = Esys[i];

  return E;
}


void HistoryPoint::extract (vector<AuxField*>& u  ,
			    real_t*              tgt) const
// ---------------------------------------------------------------------------
// Load tgt with information extracted from each AuxField in u, tgt is
// assumed to have sufficient storage to suit.
// ---------------------------------------------------------------------------
{
  register int_t      i;
  const int_t         N   = u.size();
  static const int_t  FIX = Femlib::ivalue ("LDV_FIXED");
  static const real_t w1  = Femlib::value  ("OMEGA_1");
  const real_t        wt = w1 * Femlib::value ("t");
  const real_t        wp = w1 * _y;

  if (FIX) {
    for (i = 0; i < N; i++)
      tgt[i] = u[i] -> probe (_E, _r, _s, _z - wt);
    tgt[2] += wp;
  } else
    for (i = 0; i < N; i++) tgt[i] = u[i] -> probe (_E, _r, _s, _z);
}
