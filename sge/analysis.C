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
// History points are also set up here.  They are nominated in the
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

  // -- Set up for output of plane-integral data.
  
  strcat (strcpy (str, src -> name), ".int");
  ROOTONLY {

    const real lz = Femlib::value ("TWOPI/ BETA");
    const int  nz = Geometry::nZ();
    int_strm.open (str, ios::out); 
    int_strm << "#  " << lz << "  " << nz  << "     : Lz, Nz" << endl
	     << "#  Time        Z-location         Amount"    << endl
	     << "# --------------------------------------"    << endl;
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
    const integer     NZ = Geometry::nZProc();
    const integer     NP = Geometry::nProc();
    const real        dz = Femlib::value ("TWOPI / BETA") / (NZ * NP);
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

    // -- Output plane-integral data.

    tmp.setSize (NZ*NP);

    for (i = 0; i < NZ; i++) tmp[i] = src -> u[0] -> integral (i);
  
    if (NP > 1) {
      ROOTONLY
	for (j = 1; j < NP; j++) 
	  Femlib::recv (tmp() + j * NZ, NZ, j);
      else
	Femlib::send (tmp(), NZ, 0);
    }

    ROOTONLY {
      // -- Rearrange data for, compute inverse DFT.  Output.

      const integer nz = NZ * NP;
      vector<real>  wtab (2 * nz + 15);
      
      wtab[0] = tmp[1];
      for (i = 2; i < nz - 1; i++) tmp[i - 1] = tmp[i];
      tmp[nz - 1] = wtab[0];

      Femlib::rffti (nz, wtab());
      Femlib::rfftb (nz, tmp(), wtab());

      for (i = 0; i < nz; i++)
	int_strm << setw(10) << src -> time 
		 << setw(15) << i * dz
		 << setw(15) << tmp[i]
		 << endl;

      int_strm.flush();
    }
  }

  // -- Field and statistical dumps.

  src -> dump ();
}

