///////////////////////////////////////////////////////////////////////////////
//This version of analysis.C is specialized so that it computes and
//prints out forces exerted on "wall" boundary group in non-Newtonian
//flows.  The viscosity is a function of space and time.
//
// Copyright (C) 1999 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////
 
#include <nnewt.h>

nnewtAnalyser::nnewtAnalyser (Domain* D    ,
			      FEML*   feml) : Analyser (D, feml)
// ---------------------------------------------------------------------------
// Extensions to Analyser class.
// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const char routine[] = "nnewtAnalyser::nnewtAnalyser";
    char       str[StrMax];

    // -- Open state-variable file.

    flx_strm.open (strcat (strcpy (str, _src -> name), ".flx"));
    if (!flx_strm) message (routine, "can't open flux file",  ERROR);

    flx_strm << "# nnewt state information file"      << endl;
    flx_strm << "# Step Time [Fpre Fvis Ftot]-axis" << endl;
    flx_strm << "# -------------------------------" << endl;
  }
}


void nnewtAnalyser::analyse (AuxField** work, AuxField** temp, AuxField* NNV)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  const integer DIM = Geometry::nDim();
  integer i;
  
  Analyser::analyse (work);
  
  const integer periodic = !(_src->step % (integer)Femlib::value("IO_HIS")) ||
                           !(_src->step % (integer)Femlib::value("IO_FLD"));
  const integer final    =   _src->step ==(integer)Femlib::value("N_STEP");
  const integer state    = periodic || final;
  
  if (!state) return;

  if (Geometry::system() == Geometry::Cylindrical) {
    flx_strm << "*** Cylindrical Analysis Not Yet Coded"      << endl;
    return;
  }

  Vector pfor, vfor, tfor;
  char   s[StrMax];
  
  // Need to add KINVIS (and then subtract later) because of the way
  // in which the viscous terms are handled.  Note that KINVIS has
  // been replaced by REFVIS (i.e. reference viscosity) in NS.C
  
  ROOTONLY { NNV -> addToPlane (0, Femlib::value ("KINVIS")); }
  
  NNV -> transform (INVERSE);

  //
  // First do the 2-D components
  //

  *(work[0]) = *(_src -> u[0]);
  *(temp[0]) = *(_src -> u[1]);
  *(work[1]) = *(work[0]);
  *(temp[1]) = *(temp[0]);

  for ( i=0 ; i<2 ; i++ ) {
    (*(work[i])) . gradient (i);
    (*(temp[i])) . gradient (i);
    work[i] -> transform (INVERSE);
    temp[i] -> transform (INVERSE);
    *(work[i]) *= *NNV;
    *(temp[i]) *= *NNV;
    work[i] -> transform (FORWARD);
    temp[i] -> transform (FORWARD);
  }

  ROOTONLY {

    tfor.z = pfor.z = vfor.z = 0.0;
    pfor   = Field::normTraction (_src -> u[DIM]);
    vfor   = Field::xytangTractionNN (_src -> u[0], work, temp);
    tfor.x = pfor.x + vfor.x;
    tfor.y = pfor.y + vfor.y;
    tfor.z = pfor.z;
  }

  // Now do the Z-component if needed

  if (DIM == 3) {
    
    *(work[0]) = *(_src -> u[2]);
    *(work[1]) = *(work[0]);
    
    for ( i=0 ; i<2 ; i++ ) {
      (*(work[i])) . gradient (i);
      work[i] -> transform (INVERSE);
      *(work[i]) *= *NNV;
      work[i] -> transform (FORWARD);
    }
    
    ROOTONLY {
      
      vfor.z  = (Field::ztangTractionNN (_src -> u[2], work)).z;
      tfor.z += vfor.z;

    }

  }
  
  ROOTONLY { sprintf (s,
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

  NNV -> transform (FORWARD);

  ROOTONLY { NNV -> addToPlane (0, Femlib::value ("-KINVIS")); }


}
