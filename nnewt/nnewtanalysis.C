///////////////////////////////////////////////////////////////////////////////
// This version of analysis.C is specialized so that it computes and
// prints out forces exerted on "wall" boundary group in non-Newtonian
// flows.  The viscosity is a function of space and time.
//
// Copyright (c) 1999 <--> $Date$, Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";
 
#include <nnewt.h>

nnewtAnalyser::nnewtAnalyser (Domain* D   ,
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


void nnewtAnalyser::analyse (AuxField** work,
			     AuxField** temp,
			     AuxField*  NNV )
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  const char routine[] = "nnewtAnalyser::analyse";
  const int_t DIM = Geometry::nDim();
  int_t       i;

  ROOTONLY NNV -> addToPlane (0, Femlib::value ("KINVIS"));  

  
  const bool periodic = !(_src->step %  Femlib::ivalue("IO_HIS")) ||
                        !(_src->step %  Femlib::ivalue("IO_FLD"));
  const bool final    =   _src->step == Femlib::ivalue("N_STEP");
  const bool state    = periodic || final;
  
  if (state) {

    if (Geometry::cylindrical())
      flx_strm << "Cylindrical analysis for tractions not implemented" << endl;

    // -- We are going to work out loads on walls:

    Vector pfor, vfor, tfor;
    char   s[StrMax];
  
    NNV -> transform (INVERSE);

    // -- First do the 2-D components:

    *(work[0]) = *(_src -> u[0]);
    *(temp[0]) = *(_src -> u[1]);
    *(work[1]) = *(work[0]);
    *(temp[1]) = *(temp[0]);

    for (i = 0; i < 2; i++) {
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

    // -- Now do the Z-component if needed. 

    if (DIM == 3) {
      *(work[0]) = *(_src -> u[2]);
      *(work[1]) = *(work[0]);

      for (i = 0; i < 2; i++) {
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
  
    ROOTONLY {
      sprintf (s,
	       "%6d %#10.6g "
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

  // -- This is done after all the wall computations in nnewt because for
  //    AVERAGE = 3 the pressure (needed for wall loads) will get destroyed.
  //    However, any HISTORY POINT DATA FOR PRESSURE will then also be wrong.

  Analyser::analyse (work, temp);
}
