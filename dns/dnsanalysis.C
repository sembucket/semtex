///////////////////////////////////////////////////////////////////////////////
// This version of analysis.C is specialized so that it computes and
// prints out forces exerted on "wall" boundary group.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";
 
#include "dns.h"


DNSAnalyser::DNSAnalyser (Domain* D   ,
			  BCmgr*  B   ,
			  FEML*   feml) :
// ---------------------------------------------------------------------------
// Extensions to Analyser class.
// ---------------------------------------------------------------------------
  Analyser (D, feml),
  _wss (Femlib::ivalue ("IO_WSS") && B -> nWalls())
{
  const char routine[] = "DNSAnalyser::DNSAnalyser";
  char       str[StrMax];

  ROOTONLY {
    // -- Open state-variable file.

    _flx_strm.open (strcat (strcpy (str, _src -> name), ".flx"));
    if (!_flx_strm) message (routine, "can't open flux file", ERROR);

    _flx_strm << "# DNS state information file"      << endl;
    _flx_strm << "# Step Time [Fpre Fvis Ftot]-axis" << endl;
    _flx_strm << "# -------------------------------" << endl;
  }

  if (_wss) {
    // -- Set up to compute wall shear stresses.    
    
    const int_t npr = Geometry::nProc();
    const int_t np  = Geometry::nP();
    const int_t nz  = Geometry::nZProc();

    // -- Allocate storage area: 5 for 2 normal components + 3 tangential.
    
    _nwall = B -> nWalls();
    _nline = np * _nwall;
    _npad  = 5  * _nline;
    
    // -- Round up length for Fourier transform/exchange.

    if   (npr > 1) _npad += 2 * npr - _npad % (2 * npr);
    else           _npad += _npad % 2;

    _work.resize (_npad);

    // -- Open file.

    ROOTONLY {
      _wss_strm.open (strcat (strcpy (str, _src -> name), ".wss"));
      if (!_wss_strm) message (routine, "can't open WSS file", ERROR);
    }
  }
}


void DNSAnalyser::analyse (AuxField** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "DNSAnalyser::analyse";
  const int_t DIM = Geometry::nDim();
  bool        periodic = !(_src->step %  Femlib::ivalue ("IO_HIS")) ||
                         !(_src->step %  Femlib::ivalue ("IO_FLD"));
  bool        final    =   _src->step == Femlib::ivalue ("N_STEP");
  bool        state    = periodic || final;

  Analyser::analyse (work);

  if (state) ROOTONLY {
    Vector pfor, vfor, tfor;
    char   s[StrMax];

    if (DIM == 3) {
      pfor   = Field::normTraction (_src -> u[3]);
      vfor   = Field::tangTraction (_src -> u[0], _src -> u[1], _src->u[2]);
      tfor.x = pfor.x + vfor.x;
      tfor.y = pfor.y + vfor.y;
      tfor.z = pfor.z + vfor.z;
    } else {
      pfor   = Field::normTraction (_src -> u[2]);
      vfor   = Field::tangTraction (_src -> u[0], _src -> u[1]);
      tfor.x = pfor.x + vfor.x;
      tfor.y = pfor.y + vfor.y;
      tfor.z = pfor.z = vfor.z = 0.0;
    }

    sprintf (s,
	     "%#6d %#10.6g "
	     "%#10.6g %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g "
	     "%#10.6g %#10.6g %#10.6g",
	     _src -> step, _src -> time,
	     pfor.x,   vfor.x,   tfor.x,
	     pfor.y,   vfor.y,   tfor.y,
	     pfor.z,   vfor.z,   tfor.z);

    _flx_strm << s << endl;
  }

  periodic = !(_src->step %  Femlib::ivalue ("IO_WSS")) ||
             !(_src->step %  Femlib::ivalue ("IO_FLD"));
  state    = periodic || final;

  if (_wss && state) {
    const int_t    nP  = Geometry::nP();
    const int_t    nZ  = Geometry::nZ();
    const int_t    nZP = Geometry::nZProc();
    const int_t    nPR = Geometry::nProc();
    const int_t    nPP = _npad / nPR;
    int_t          i, j, k;
    real_t*        plane;
    vector<real_t> buffer (_nline);

    // -- Load the local storage area.

    Veclib::zero (_work.size(), &_work[0], 1);

    if (DIM == 3 || _src -> nField() == 4)
      Field::traction (&_work[0], &_work[_nline], &_work[2*_nline],
		       &_work[3*_nline], &_work[4*_nline], _nwall,
		       _src -> u[3], _src -> u[0], _src -> u[1], _src -> u[2]);
    else
      Field::traction (&_work[0], &_work[_nline], &_work[2*_nline],
		       &_work[3*_nline], &_work[4*_nline], _nwall,
		       _src -> u[2], _src -> u[0], _src -> u[1]);

    // -- Inverse Fourier transform (like Field::bTransform).

    if (nPR == 1) {
      if (nZ > 1)
	if (nZ == 2)
	  Veclib::copy (_npad, &_work[0], 1, &_work[_npad], 1);
	else
	  Femlib::DFTr (&_work[0], nZ, _npad, INVERSE);
    } else {
      Femlib::exchange (&_work[0], nZP, nP,  FORWARD);
      Femlib::DFTr     (&_work[0], nZ,  nPP, INVERSE);
      Femlib::exchange (&_work[0], nZP, nP,  INVERSE);
    }

    // -- Write to file.

    // -- Header: this will be a lot like a standard header.

    ROOTONLY {
      const char *hdr_fmt[] = { 
	"%-25s "                "Session\n",
	"%-25s "                "Created\n",
	"%-5d %-5d %-5d      "  "Np, Nz, Elements\n",
	"%-25d "                "Step\n",
	"%-25.6g "              "Time\n",
	"%-25.6g "              "Time step\n",
	"%-25.6g "              "Kinvis\n",
	"%-25.6g "              "Beta\n",
	"%-25s "                "Fields written\n",
	"%-25s "                "Format\n"
      };

      char   s1[StrMax], s2[StrMax];
      time_t tp (time (0));

      strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));

      sprintf (s1, hdr_fmt[0], _src->name);               _wss_strm << s1;
      sprintf (s1, hdr_fmt[1], s2);                       _wss_strm << s1;
      sprintf (s1, hdr_fmt[2], nP, nZ, _nwall);           _wss_strm << s1;
      sprintf (s1, hdr_fmt[3], _src->step);               _wss_strm << s1;
      sprintf (s1, hdr_fmt[4], _src->time);               _wss_strm << s1;
      sprintf (s1, hdr_fmt[5], Femlib::value ("D_T"));    _wss_strm << s1;
      sprintf (s1, hdr_fmt[6], Femlib::value ("KINVIS")); _wss_strm << s1;
      sprintf (s1, hdr_fmt[7], Femlib::value ("BETA"));   _wss_strm << s1;
      sprintf (s1, hdr_fmt[8], "NxNyTxTyTz");             _wss_strm << s1;
      sprintf (s2, "binary "); Veclib::describeFormat (s2 + strlen (s2));
      sprintf (s1, hdr_fmt[9], s2);                       _wss_strm << s1;

      if (!_wss_strm) message (routine, "failed writing WSS header", ERROR);
      _wss_strm << flush;
      
    }

    // -- Data.

    if (nPR > 1) {		// -- Parallel.
      for (j = 0; j < 5; j++)	// -- Reminder: there are 5 components.
	ROOTONLY {
	  for (i = 0; i < nZP; i++) {
	    plane = &_work[i*_npad + j*_nline];
	    _wss_strm.write(reinterpret_cast<char*>(plane),
			    static_cast<int_t>(_nline * sizeof (real_t))); 
	    if (_wss_strm.bad())
	      message (routine, "unable to write binary output", ERROR);
	  }
	  for (k = 1; k < nPR; k++)
	    for (i = 0; i < nZP; i++) {
	      Femlib::recv (&buffer[0], _nline, k);
	      _wss_strm.write(reinterpret_cast<char*>(&buffer[0]),
			      static_cast<int_t>(_nline * sizeof (real_t))); 
	      if (_wss_strm.bad()) 
		message (routine, "unable to write binary output", ERROR);
	    }
        } else			// -- Not on root process.
	  for (i = 0; i < nZP; i++) {
	    plane = &_work[i*_npad + j*_nline];
	    Femlib::send (plane, _nline, 0);
	  }
    } else			// -- Serial.
      for (j = 0; j < 5; j++)
	for (i = 0; i < nZ; i++) {
	  plane = &_work[i*_npad + j*_nline];
	  _wss_strm.write (reinterpret_cast<char*>(plane),
			   static_cast<int_t>(_nline * sizeof (real_t))); 
	  if (_wss_strm.bad())
	    message (routine, "unable to write binary output", ERROR);
	}
  }
}
