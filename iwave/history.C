///////////////////////////////////////////////////////////////////////////////
// history.C
//
// Copyright (C) 1994, 1999 Hugh Blackburn
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

#include <Sem.h>


const Element* HistoryPoint::locate (const real        x   ,
				     const real        y   ,
				     vector<Element*>& Esys,
				     real&             r   ,
				     real&             s   )
// ---------------------------------------------------------------------------
// Static class member function which tries to locate x, y, point within
// an element E, and return its location in r, s coordinates within E.
// ---------------------------------------------------------------------------
{
  register integer i;
  const integer    NEL    = Esys.getSize();
  const integer    guess = 1;
  const Element*   E;

  E = 0;
  for (i = 0; i < NEL; i++) {
    r = s = 0.0;
    if (Esys[i] -> locate (x, y, r, s, guess)) {
      E = Esys[i];
      break;
    }
  }

  return E;
}


void HistoryPoint::extract (vector<AuxField*>& u  ,
			    real*              tgt) const
// ---------------------------------------------------------------------------
// Load tgt with information extracted from each AuxField in u, tgt is
// assumed to have sufficient storage to suit.
// ---------------------------------------------------------------------------
{
  register integer i;
  const integer    N = u.getSize();

  static const integer FIX = (integer) Femlib::value ("LDV_FIXED");
  static const real    w1  =           Femlib::value ("OMEGA_1");

  const real wt = w1 * Femlib::value ("t");
  const real wp = w1 * _y;

  if (FIX) {
    for (i = 0; i < N; i++)
      tgt[i] = u[i] -> probe (_E, _r, _s, _z - wt);
    tgt[2] += wp;
  } else
    for (i = 0; i < N; i++) tgt[i] = u[i] -> probe (_E, _r, _s, _z);
}

