///////////////////////////////////////////////////////////////////////////////
// analysis.C: implement Analyser class for NS-type solvers.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include "Sem.h"


Analyser::Analyser (Domain* D   ,
		    FEML*   file) :
// ---------------------------------------------------------------------------
// Set up.
//
// History points are also set up with.  They are nominated in the
// optional HISTORY section of the session file.  Output is to
// session.his.
// ---------------------------------------------------------------------------
  src (D)
{
  const char routine[] = "Analyser::Analyser";
  char       str[StrMax];

  cout << setprecision (6);

  // -- Set up for history points: open files, create points.

  if (file -> seek ("HISTORY")) {
    integer              i, id, num = 0;
    const integer        NH = file -> attribute ("HISTORY", "NUMBER");
    const Element*       E;
    HistoryPoint*        H;
    Stack<HistoryPoint*> stack;
    real                 r, s, x, y, z;
    
    for (i = 0; i < NH; i++) {
      file -> stream() >> id >> x >> y >> z;
      if (E = HistoryPoint::locate (x, y, D -> elmt, r, s)) {
	H = new HistoryPoint (id, E, r, s, z);
	stack.push (H);
	num++;
      } else {
	sprintf (str, "History point at (%f, %f, %f) not in mesh", x, y, z);
	message (routine, str, WARNING);
      }
    }
    
    history.setSize (num);
    while (num--) history[num] = stack.pop();

    ROOTONLY {
      his_strm.open (strcat (strcpy (str, src -> name), ".his"));
      his_strm.setf (ios::scientific, ios::floatfield);
      his_strm.precision (6);
      if (!his_strm) message (routine, "can't open history file", ERROR);
    }
  }
}


void Analyser::analyse ()
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  const integer verbose = (integer) Femlib::value ("VERBOSE");
  const integer cflstep = (integer) Femlib::value ("IO_CFL");

  // -- Step-by-step updates.

  ROOTONLY cout << "Step: " << src->step << "  Time: " << src->time << endl;

  // -- Periodic dumps and global information.
  
  const integer periodic = !(src->step %  (integer) Femlib::value("IO_HIS")) ||
                           !(src->step %  (integer) Femlib::value("IO_FLD"));
  const integer final    =   src->step == (integer) Femlib::value("N_STEP");
  const integer state    = periodic || final;

  if (state) {

    // -- Output history point data.
      
    register integer  i, j;
    const integer     NH = history.getSize();
    const integer     NF = src -> u.getSize();
    HistoryPoint*     H;
    vector<real>      tmp (NF);
    vector<AuxField*> u   (NF);

    for (i = 0; i < NF; i++)
      u[i] = src -> u[i];

    for (i = 0; i < NH; i++) {
      H = history[i];

      H -> extract (u, tmp());

      ROOTONLY {
	his_strm << setw(4) << H->ID() << " " << setw(14) << src->time << " ";
	for (j = 0; j < NF; j++) his_strm << setw(15) << tmp[j];
	his_strm << endl;
      }
    }
  }

  // -- Field and statistical dumps.

  src -> dump ();
}

