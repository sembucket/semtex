///////////////////////////////////////////////////////////////////////////////
// This version of analysis.C is specialized so that it computes and
// prints out the temperature flux on "wall" boundary group.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////
 
#include "scat.h"


ScatAnalyser::ScatAnalyser (Domain* D   ,
			    FEML*   feml) :
// ---------------------------------------------------------------------------
// Extensions to Analyser class.
// ---------------------------------------------------------------------------
  Analyser (D, feml)
{
  ROOTONLY {
    const char routine[] = "BuoyAnalyser::BuoyAnalyser";
    char       str[StrMax];

    // -- Open state-variable file.

    flx_strm.open (strcat (strcpy (str, src -> name), ".flx"));
    if (!flx_strm) message (routine, "can't open flux file",  ERROR);

    flx_strm << "# Buoy state information file"          << endl;
    flx_strm << "# Step Time Flux [Fpre Fvis Ftot]-axis" << endl;
    flx_strm << "# ------------------------------------" << endl;
  }
}


void ScatAnalyser::analyse (AuxField** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  const integer DIM = Geometry::nDim();

  Analyser::analyse (work);

  ROOTONLY {
    const integer periodic = !(src->step% (integer)Femlib::value ("IO_HIS")) ||
                             !(src->step% (integer)Femlib::value ("IO_FLD"));
    const integer final    =   src->step==(integer)Femlib::value ("N_STEP");
    const integer state    = periodic || final;

    if (!state) return;

    real   flux;
    Vector pfor, vfor, tfor;
    char   s[StrMax];

    flux   = Field::flux            (src -> u[DIM]);
    pfor   = Field::normalTraction  (src -> u[DIM + 1]);
    vfor   = Field::tangentTraction (src -> u[0], src-> u[1]);
    tfor.x = pfor.x + vfor.x;
    tfor.y = pfor.y + vfor.y;

    sprintf (s,
	     "%#6d %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g %#10.6g %#10.6g %#10.6g",
	     src -> step, src -> time, flux  ,
	     pfor.x,   vfor.x,   tfor.x,
	     pfor.y,   vfor.y,   tfor.y);

    flx_strm << s << endl;
  }
}
