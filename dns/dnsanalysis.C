///////////////////////////////////////////////////////////////////////////////
// This version of analysis.C is specialized so that it computes and
// prints out forces exerted on "wall" boundary group.
//
// Copyright (C) 1994, 2000 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////
 
#include <dns.h>


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

    flx_strm.open (strcat (strcpy (str, src -> name), ".flx"));
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
    const integer periodic = !(src->step % (integer)Femlib::value("IO_HIS")) ||
                             !(src->step % (integer)Femlib::value("IO_FLD"));
    const integer final    =   src->step ==(integer)Femlib::value("N_STEP");
    const integer state    = periodic || final;

    if (!state) return;

    Vector pfor, vfor, tfor;
    char   s[StrMax];

    if (DIM == 3) {
      pfor   = Field::normalTraction  (src -> u[3]);
      vfor   = Field::tangentTraction (src -> u[0], src -> u[1], src -> u[2]);
      tfor.x = pfor.x + vfor.x;
      tfor.y = pfor.y + vfor.y;
      tfor.z = pfor.z + vfor.z;
    } else {
      pfor   = Field::normalTraction  (src -> u[2]);
      vfor   = Field::tangentTraction (src -> u[0], src -> u[1]);
      tfor.x = pfor.x + vfor.x;
      tfor.y = pfor.y + vfor.y;
      tfor.z = pfor.z = vfor.z = 0.0;
    }

    sprintf (s,
	     "%#6d %#10.6g "
	     "%#10.6g %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g",
	     src -> step, src -> time,
	     pfor.x,   vfor.x,   tfor.x,
	     pfor.y,   vfor.y,   tfor.y,
	     pfor.z,   vfor.z,   tfor.z);

    flx_strm << s << endl;
  }
}
