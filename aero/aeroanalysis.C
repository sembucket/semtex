///////////////////////////////////////////////////////////////////////////////
// aeroanalysis.C: implement AeroAnalyser class,
// an extension to Analyser class.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "aero.h"


AeroAnalyser::AeroAnalyser (Domain* D   ,
			    FEML*   feml,
			    Body*   B   ) : Analyser (D, feml), body (B)
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

    sta_strm.open (strcat (strcpy (str, src -> name), ".sta"));

    if (!sta_strm) message (routine, "can't open state file", ERROR);

    sta_strm << "# Aero state information file"                 << endl;
    sta_strm << "# Step Time [pos vel acc Fpre Fvis Ftot]-axis" << endl;
    sta_strm << "# -------------------------------------------" << endl;

#if defined(FORCES)
    // -- Open force file.

    for_strm.open (strcat (strcpy (str, src -> name), ".for"));

    if (!for_strm) message (routine, "can't open force file", ERROR);

    for_strm << "# Aero force information file" << endl;
    for_strm << "# Time ZPos  Fx  Fy  Fz"       << endl;
    for_strm << "# ---------------------------" << endl;
#endif
  }
}


void AeroAnalyser::analyse (AuxField** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.
// ---------------------------------------------------------------------------
{
  Analyser::analyse (work);

  ROOTONLY {
    const integer periodic = !(src->step% (integer)Femlib::value ("IO_HIS")) ||
                             !(src->step% (integer)Femlib::value ("IO_FLD"));
    const integer final    =   src->step==(integer)Femlib::value ("N_STEP");
    const integer state    = periodic || final;

    if (!state) return;

    sta_strm << src -> step << " " << src -> time << body << endl;
  }

#if defined(FORCES)
  const integer forstep = (integer) Femlib::value ("IO_HIS");

  if (forstep && !(src->step % forstep) && Geometry::nDim() == 3) forceDist();
#endif 
}


void AeroAnalyser::forceDist ()
// ---------------------------------------------------------------------------
// Compute forces at each axial location, output on root processor.
// Designed only to be called for 3D.
// ---------------------------------------------------------------------------
{
  register integer i, k;
  const integer    nZ     = Geometry::nZ();
  const integer    nProc  = Geometry::nProc();
  const integer    nZProc = Geometry::nZProc();
  Field*           master = src -> u[0];
  vector<real>     work (5 * nZProc);
  register real    *px = &work[0],
                   *py = px + nZProc,
                   *vx = py + nZProc,
                   *vy = vx + nZProc,
                   *vz = vy + nZProc;

  // -- Fill local copy of force data.

  Veclib::zero (5 * nZProc, px, 1);
  master -> normTractionV (px, py,     src -> u[3]);
  master -> tangTractionV (vx, vy, vz, src -> u[0], src -> u[1], src -> u[2]);
  Veclib::vadd (nZProc, vx, 1, px, 1, vx, 1);
  Veclib::vadd (nZProc, vy, 1, py, 1, vy, 1);

  // -- Order data, inverse Fourier transform, print up.

  if (nProc > 1) {

    ROOTONLY {
      vector<real>  force (4 * nZ);	// -- Even number forced by DFTr.
      register real *x = &force[0], *y = x + 1, *z = y + 1, *t = z + 1;
      const real    dz = Femlib::value ("TWOPI / BETA / N_Z");

      for (i = 0; i < nZProc; i++) {
	x[4 * i] = vx[i];
	y[4 * i] = vy[i];
	z[4 * i] = vz[i];
	t[4 * i] = 0.0;
      }
	
      for (k = 1; k < nProc; k++) {
	Femlib::recv (vx, 3 * nZProc, k);
	for (i = 0; i < nZProc; i++) {
	  x[4 * (i + k * nZProc)] = vx[i];
	  y[4 * (i + k * nZProc)] = vy[i];
	  z[4 * (i + k * nZProc)] = vz[i];
	  t[4 * (i + k * nZProc)] = 0.0;
	}
      }
      
      Femlib::DFTr (x, nZ, 4, -1);

      for (k = 0; k < nZ; k++)
	for_strm << setw(10) << src -> time 
		 << setw(15) << k * dz
		 << setw(15) << x[4 * k]
		 << setw(15) << y[4 * k]
		 << setw(15) << z[4 * k]
		 << endl;

      for_strm.flush();

    } else
      Femlib::send (vx, 3 * nZProc, 0);

  } else {			// -- Serial version.

    vector<real>  force (4 * nZ);
    register real *x = &force[0], *y = x + 1, *z = y + 1, *t = z + 1;
    const real    dz = Femlib::value ("TWOPI / BETA / N_Z");

    for (k = 0; k < nZ; k++) {
      x[4 * k] = vx[k];
      y[4 * k] = vy[k];
      z[4 * k] = vz[k];
      t[4 * k] = 0.0;
    }
	
    Femlib::DFTr (x, nZ, 4, -1);

    for (k = 0; k < nZ; k++)
      for_strm << setw(10) << src -> time 
	       << setw(15) << k * dz
	       << setw(15) << x[4 * k]
	       << setw(15) << y[4 * k]
	       << setw(15) << z[4 * k]
	       << endl;

    for_strm.flush();
  }
}
