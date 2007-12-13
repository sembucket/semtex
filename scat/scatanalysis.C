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
    const char routine[] = "ScatAnalyser::ScatAnalyser";
    char       str[StrMax];

    // -- Open state-variable file.

    flx_strm.open (strcat (strcpy (str, _src -> name), ".flx"));
    if (!flx_strm) message (routine, "can't open flux file",  ERROR);

    flx_strm << "# Scat state information file"          << endl;
    flx_strm << "# Step Time Flux [Fpre Fvis Ftot]-axis" << endl;
    flx_strm << "# ------------------------------------" << endl;
  }
}


void ScatAnalyser::analyse (AuxField** wrk1, AuxField** wrk2)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  const int_t DIM = Geometry::nDim();

  Analyser::analyse (wrk1, wrk2);

  ROOTONLY {
    const bool periodic = !(_src->step  % Femlib::ivalue ("IO_HIS")) ||
                          !(_src->step  % Femlib::ivalue ("IO_FLD"));
    const bool final    =   _src->step == Femlib::ivalue ("N_STEP");
    const bool state    = periodic || final;

    if (!state) return;

    real_t flux;
    Vector pfor, vfor, tfor;
    char   s[StrMax];

    flux   = Field::scalarFlux   (_src -> u[DIM]);
    pfor   = Field::normTraction (_src -> u[DIM + 1]);
    vfor   = Field::tangTraction (_src -> u[0], _src-> u[1]);
    tfor.x = pfor.x + vfor.x;
    tfor.y = pfor.y + vfor.y;

    sprintf (s,
	     "%#6d %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g %#10.6g %#10.6g %#10.6g",
	     _src -> step, _src -> time, flux  ,
	     pfor.x,   vfor.x,   tfor.x,
	     pfor.y,   vfor.y,   tfor.y);

    flx_strm << s << endl;
  }
}
