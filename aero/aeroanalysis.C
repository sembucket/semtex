///////////////////////////////////////////////////////////////////////////////
// aeroanalysis.C: implement AeroAnalyser class,
// an extension to Analyser class.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <aero.h>


AeroAnalyser::AeroAnalyser (Domain& D   ,
			    FEML&   feml,
			    Body&   B   ) : Analyser (D, feml), body (B)
// ---------------------------------------------------------------------------
// Extensions to Analyser class.
//
// NB: at present particles are spawned and move in relation to
// moving reference frame.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AeroAnalyser::AeroAnalyser";
  char       str[StrMax];

  // -- Open state-variable file.

  sta_strm.open (strcat (strcpy (str, src.name), ".sta"));

  if (!sta_strm) message (routine, "can't open state file", ERROR);

  sta_strm << "# Aero state information file"                 << endl;
  sta_strm << "# Step Time [pos vel acc Fpre Fvis Ftot]-axis" << endl;
  sta_strm << "# -------------------------------------------" << endl;
}


void AeroAnalyser::analyse (AuxField*** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  Analyser::analyse (work);

  sta_strm << src.step << " " << src.time << body << endl;
}
