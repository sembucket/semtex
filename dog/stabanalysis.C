///////////////////////////////////////////////////////////////////////////////
// This version of analysis.C is specialized so that it
// prints out base history point information.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
// RCSid[] = "$Id$";
//
///////////////////////////////////////////////////////////////////////////////
 
#include <stab.h>


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

  if (file -> seek ("BASE")) {
    integer              i, id, num = 0;
    const integer        NBH = file -> attribute ("BASE_HIST", "NUMBER");
    const Element*       EB;
    HistoryPoint*        HB;
    Stack<HistoryPoint*> base_stack;
    real                 r, s, x, y, z;
    
    for (i = 0; i < NBH; i++) {
      file -> stream() >> id >> x >> y >> z;
      if ((EB = HistoryPoint::locate (x, y, D -> elmt, r, s))) {
	HB = new HistoryPoint (id, EB, r, s, z);
	base_stack.push (HB);
	num++;
      } else {
	sprintf (str, "History point at (%f, %f, %f) not in mesh", x, y, z);
	message (routine, str, WARNING);
      }
    }
    
    base_history.setSize (num);
    while (num--) base_history[num] = base_stack.pop();
    
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

  const integer periodic =
    !(src->step %  static_cast<integer>(Femlib::value("IO_HIS"))) ||
    !(src->step %  static_cast<integer>(Femlib::value("IO_FLD")));
  const integer final    =
    src->step == static_cast<integer>(Femlib::value("N_STEP"));
  const integer state    =
    periodic || final;

  if (!state) return;

  // -- Output BASE history point data.

  register integer  j, k;      
  const integer     NBH = base_history.getSize();
  const integer     NBF = 2;  // number of base fields = 2 (UV)
  HistoryPoint*     HB;
  vector<real>      tmp_B (NBF);
  vector<AuxField*> U     (NBF);
  
  for (k = 0; k < NBF; k++)
    U[k] = src -> U[k];

  for (k = 0; k < NBH; k++) {
    HB = base_history[k];
    
    HB -> extract (U, tmp_B());

    bhs_strm << setw(4) << HB->ID() << " " << setw(14) << src->time << " ";
    for (j = 0; j < NBF; j++) bhs_strm << setw(15) << tmp_B[j];
    bhs_strm << endl;
  }
}
