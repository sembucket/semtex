///////////////////////////////////////////////////////////////////////////////
// This version of analysis.C is specialized so that it computes and
// prints out the temperature flux on "wall" boundary group.
///////////////////////////////////////////////////////////////////////////////
 
static char
RCSid[] = "$Id$";

#include <buoy.h>


BuoyAnalyser::BuoyAnalyser (Domain& D   ,
			    FEML&   feml) : Analyser (D, feml)
// ---------------------------------------------------------------------------
// Extensions to Analyser class.
// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const char routine[] = "BuoyAnalyser::BuoyAnalyser";
    char       str[StrMax];

    // -- Open state-variable file.

    flx_strm.open (strcat (strcpy (str, src.name), ".flx"));

    if (!flx_strm) message (routine, "can't open flux file", ERROR);
  }
}


void BuoyAnalyser::analyse (AuxField*** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
 const integer DIM = Geometry::nDim();

 Analyser::analyse (work);

  ROOTONLY {
    const integer periodic = !(src.step % (integer)Femlib::value ("IO_HIS")) ||
                             !(src.step % (integer)Femlib::value ("IO_FLD"));
    const integer final    =   src.step ==(integer)Femlib::value ("N_STEP");
    const integer state    = periodic || final;

    if (!state) return;

    flx_strm << src.step << " " << src.time;
    flx_strm << Field::flux(src.u[DIM+1]) << endl;
  }
}
