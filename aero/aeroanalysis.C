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
  ROOTONLY {
    const char routine[] = "AeroAnalyser::AeroAnalyser";
    char       str[StrMax];

    // -- Open state-variable file.

    sta_strm.open (strcat (strcpy (str, src.name), ".sta"));

    if (!sta_strm) message (routine, "can't open state file", ERROR);

    sta_strm << "# Aero state information file"                 << endl;
    sta_strm << "# Step Time [pos vel acc Fpre Fvis Ftot]-axis" << endl;
    sta_strm << "# -------------------------------------------" << endl;
  }
}


void AeroAnalyser::analyse (AuxField*** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  Analyser::analyse (work);

  ROOTONLY {
    const integer periodic = !(src.step % (integer)Femlib::value ("IO_HIS")) ||
                             !(src.step % (integer)Femlib::value ("IO_FLD"));
    const integer final    =   src.step ==(integer)Femlib::value ("N_STEP");
    const integer state    = periodic || final;

    if (!state) return;

    sta_strm << src.step << " " << src.time << body << endl;
  }
}
