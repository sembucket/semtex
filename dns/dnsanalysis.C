///////////////////////////////////////////////////////////////////////////////
// This version of analysis.C is specialized so that it computes and
// prints out forces exerted on "wall" boundary group.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";
 
#include "dns.h"


DNSAnalyser::DNSAnalyser (Domain* D   ,
			  FEML*   feml) :
// ---------------------------------------------------------------------------
// Extensions to Analyser class.
// ---------------------------------------------------------------------------
  Analyser (D, feml)
{
  ROOTONLY {
    const char routine[] = "DNSAnalyser::DNSAnalyser";
    char       str[StrMax];

    // -- Open state-variable file.

    flx_strm.open (strcat (strcpy (str, _src -> name), ".flx"));
    if (!flx_strm) message (routine, "can't open flux file",  ERROR);

    flx_strm << "# DNS state information file"      << endl;
    flx_strm << "# Step Time [Fpre Fvis Ftot]-axis" << endl;
    flx_strm << "# -------------------------------" << endl;
  }
}


void DNSAnalyser::analyse (AuxField** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  const integer DIM = Geometry::nDim();

  Analyser::analyse (work);

  ROOTONLY {
    const bool periodic = !(_src->step %  Femlib::ivalue ("IO_HIS")) ||
                          !(_src->step %  Femlib::ivalue ("IO_FLD"));
    const bool final    =   _src->step == Femlib::ivalue ("N_STEP");
    const bool state    = periodic || final;

    if (!state) return;

    Vector pfor, vfor, tfor;
    char   s[StrMax];

    if (DIM == 3) {
      pfor   = Field::normalTraction  (_src -> u[3]);
      vfor   = Field::tangentTraction (_src -> u[0], _src -> u[1], _src -> u[2]);
      tfor.x = pfor.x + vfor.x;
      tfor.y = pfor.y + vfor.y;
      tfor.z = pfor.z + vfor.z;
    } else {
      pfor   = Field::normalTraction  (_src -> u[2]);
      vfor   = Field::tangentTraction (_src -> u[0], _src -> u[1]);
      tfor.x = pfor.x + vfor.x;
      tfor.y = pfor.y + vfor.y;
      tfor.z = pfor.z = vfor.z = 0.0;
    }

    sprintf (s,
	     "%#6d %#10.6g "
	     "%#10.6g %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g",
	     _src -> step, _src -> time,
	     pfor.x,   vfor.x,   tfor.x,
	     pfor.y,   vfor.y,   tfor.y,
	     pfor.z,   vfor.z,   tfor.z);

    flx_strm << s << endl;
  }
}
