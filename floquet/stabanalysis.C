///////////////////////////////////////////////////////////////////////////////
// This version of analysis.C is specialized so that it
// prints out base history point information.
//
// Copyright (C) 1994,2003 Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";
 
#include "stab.h"


StabAnalyser::StabAnalyser (Domain* D   ,
			    FEML*   file) :
// ---------------------------------------------------------------------------
// Extensions to Analyser class.
// ---------------------------------------------------------------------------
  Analyser (D, file)
{
  // -- Open base history file.
  
  const char routine[] = "StabAnalyser::StabAnalyser";
  char       str[StrMax];

  if (file -> seek ("BASE_HIST")) {
    int            i, id, num = 0;
    const int      NBH = file -> attribute ("BASE_HIST", "NUMBER");
    const Element* EB;
    HistoryPoint*  HB;
    real           r, s, x, y, z;
    
    for (i = 0; i < NBH; i++) {
      file -> stream() >> id >> x >> y >> z;
      if ((EB = HistoryPoint::locate (x, y, D -> elmt, r, s))) {
	HB = new HistoryPoint (id, EB, r, s, z);
	base_history.insert (base_history.end(), HB);
	num++;
      } else {
	sprintf (str, "History point at (%f, %f, %f) not in mesh", x, y, z);
	message (routine, str, WARNING);
      }
    }
    
    bhs_strm.open (strcat (strcpy (str, src -> name), ".bhs"));
    bhs_strm.setf (ios::scientific, ios::floatfield);
    bhs_strm.precision (6);
    if (!bhs_strm) message (routine, "can't open history file", ERROR);
  }


}


void StabAnalyser::analyse (AuxField** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  Analyser::analyse (work);

  const int periodic =
    !(src->step %  static_cast<int>(Femlib::value("IO_HIS"))) ||
    !(src->step %  static_cast<int>(Femlib::value("IO_FLD")));
  const int final    =
    src->step == static_cast<int>(Femlib::value("N_STEP"));
  const int state    =
    periodic || final;

  if (!state) return;

  // -- Output BASE history point data.

  register int      j, k;      
  const int         NBH = base_history.size();
  const int         NBF = 2;  // number of base fields = 2 (UV)
  HistoryPoint*     HB;
  vector<real>      tmp_B (NBF);
  vector<AuxField*> U     (NBF);
  
  for (k = 0; k < NBF; k++)
    U[k] = src -> U[k];

  for (k = 0; k < NBH; k++) {
    HB = base_history[k];
    
    HB -> extract (U, &tmp_B[0]);

    bhs_strm << setw(4) << HB->ID() << " " << setw(14) << src->time << " ";
    for (j = 0; j < NBF; j++) bhs_strm << setw(15) << tmp_B[j];
    bhs_strm << endl;
  }
}
