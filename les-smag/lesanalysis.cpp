///////////////////////////////////////////////////////////////////////////////
// This version of analysis.C is specialized so that it computes and
// prints out forces exerted on "wall" boundary group.
//
// Copyright (c) 1999 <--> $Date$, Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";
 
#include "les.h"


LESAnalyser::LESAnalyser (Domain* D   ,
			  FEML*   feml) :
// ---------------------------------------------------------------------------
// Extensions to Analyser class.
// ---------------------------------------------------------------------------
  Analyser (D, feml)
{
  ROOTONLY {
    const char routine[] = "LESAnalyser::LESAnalyser";
    char       str[StrMax];

    // -- Open state-variable file.

    _flx_strm.open (strcat (strcpy (str, _src -> name), ".flx"));
    if (!_flx_strm) message (routine, "can't open flux file",  ERROR);

    _flx_strm << "# LES state information file"      << endl;
    _flx_strm << "# Step Time [Fpre Fvis Ftot]-axis" << endl;
    _flx_strm << "# -------------------------------" << endl;
  }
}


void LESAnalyser::analyse (AuxField** work0,
                           AuxField** work1)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  const int_t DIM = Geometry::nDim();

  Analyser::analyse (work0, work1);

  ROOTONLY {
    const bool periodic = !(_src->step  % Femlib::ivalue ("IO_HIS")) ||
                          !(_src->step  % Femlib::ivalue ("IO_FLD"));
    const bool final    =   _src->step == Femlib::ivalue ("N_STEP");
    const bool state    = periodic || final;

    if (!state) return;

    Vector pfor, vfor, tfor;
    char   s[StrMax];

    if (DIM == 3) {
      pfor   = Field::normTraction (_src->u[3]);
      vfor   = Field::tangTraction (_src->u[0], _src->u[1], _src->u[2]);
      tfor.x = pfor.x + vfor.x;
      tfor.y = pfor.y + vfor.y;
      tfor.z = pfor.z + vfor.z;
    } else {
      pfor   = Field::normTraction (_src->u[2]);
      vfor   = Field::tangTraction (_src->u[0], _src->u[1]);
      tfor.x = pfor.x + vfor.x;
      tfor.y = pfor.y + vfor.y;
      tfor.z = pfor.z = vfor.z = 0.0;
    }

    sprintf (s,
	     "%6d %#10.6g "
	     "%#10.6g %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g",
	     _src -> step, _src -> time,
	     pfor.x,   vfor.x,   tfor.x,
	     pfor.y,   vfor.y,   tfor.y,
	     pfor.z,   vfor.z,   tfor.z);

    _flx_strm << s << endl;
  }
}