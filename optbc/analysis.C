//////////////////////////////////////////////////////////////////////////////
// analysis.C: implement Analyser class for NS-type solvers.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// This deals with output of runtime information such as step numbers,
// CFL estimation, modal energies, etc.  If set, also output history
// point and particle track information.
//
// It is assumed that the first 2 or 3 (for 3D) entries in the Domain
// u vector are velocity fields.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <unistd.h>
#include <sem.h>


Analyser::Analyser (Domain* D   ,
		    FEML*   file) :
// ---------------------------------------------------------------------------
// Set up.
//
// Fluid particle files (session.par) are dealt with here.
// Each line is of form
// #     tag  time  x     y      z
//       1    0.0   1.0   10.0   0.5.
// Output is of the same form, called session.trk.
//
// NB: Particle tracking is broken for multiprocessor application.
//
// History points are also set up here.  They are nominated in the
// optional HISTORY section of the session file.  Output is to
// session.his.
// ---------------------------------------------------------------------------
  _src (D)
{
  const char routine[] = "Analyser::Analyser";
  char       str[StrMax];
  time_t     tp (time (0));

  cout << setprecision (6);

  // -- Set up for history points: open files, create points.

  if (file -> seek ("HISTORY")) {
    int_t          i, id, num = 0;
    const int_t    NH = file -> attribute ("HISTORY", "NUMBER");
    const Element* E;
    HistoryPoint*  H;
    real_t         r, s, x, y, z;
    
    for (i = 0; i < NH; i++) {
      file -> stream() >> id >> x >> y >> z;
      if ((E = HistoryPoint::locate (x, y, D -> elmt, r, s))) {
	H = new HistoryPoint (id, E, r, s, x, y, z);
	_history.insert (_history.end(), H);
	num++;
      } else {
	sprintf (str, "History point at (%f, %f, %f) not in mesh", x, y, z);
	message (routine, str, WARNING);
      }
    }
    
    _his_strm.open (strcat (strcpy (str, _src -> name), ".his"));
    _his_strm.setf (ios::scientific, ios::floatfield);
    _his_strm.precision (6);
    if (!_his_strm) message (routine, "can't open history file", ERROR);
  }

  // -- Set up for output of modal energies every IO_CFL steps.

  _mdl_strm.open (strcat (strcpy (str, _src -> name), ".mdl"), ios::out); 
  _mdl_strm << "#     Time Mode         Energy Energy_fo_ axial_waves" << endl
	   << "# ----------------------------" << endl;

  // -- Dump run information to file.

  ofstream runfile (strcat (strcpy (str, _src -> name), ".run"), ios::out);
  gethostname (str, StrMax);
  runfile << "-- Host                    : " << str << endl;
  runfile << "   PID                     : " << getpid() << endl;

  strftime (str, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  runfile << "   Date                    : " << str << endl;

  D -> report (runfile);
  
  runfile.close();
}


void Analyser::analyse (AuxField** work,
			AuxField** not_used)
// ---------------------------------------------------------------------------
// Step-by-step processing.  If SPAWN was set, add more particles at
// original absolute positions.
// ---------------------------------------------------------------------------
{
  const int_t cflstep = Femlib::ivalue ("IO_CFL");

  // -- Run information update.
cout << "Step: " << _src -> step << "  Time: " << _src -> time << endl; 


  // -- CFL, energy, divergence information.

  if (cflstep && !(_src -> step % cflstep)) {
    estimateCFL ();
    divergence  (work);
  }
if(!(_src -> step % Femlib::ivalue ("IO_MDL")))  
 if (Femlib::ivalue ("quadrilaturalbatchelor")==1) modalEnergyquad ();
 else modalEnergy();
  // -- Periodic dumps and global information.
  
  const bool periodic = !(_src -> step %  Femlib::ivalue("IO_HIS")) ||
                        !(_src -> step %  Femlib::ivalue("IO_FLD")) ;
  const bool final    =   _src -> step == Femlib::ivalue("N_STEP");
  const bool state    = periodic || final;

  if (state) {
    // -- Output history point data.
      
    register int_t    i, j;
    const int_t       NH = _history.size();
    const int_t       NF = _src -> u.size();
    HistoryPoint*     H;
    vector<real_t>    tmp (NF);
    vector<AuxField*> u   (NF);

    for (i = 0; i < NF; i++)
      u[i] = _src -> u[i];

    for (i = 0; i < NH; i++) {
      H = _history[i];

      H -> extract (u, &tmp[0]);

      _his_strm << setw(4) << H->ID() << " " << setw(14) << _src->time << " ";
      for (j = 0; j < NF; j++) _his_strm << setw(15) << tmp[j];
      _his_strm << endl;
    }
  }

  _src -> dump();
}

void Analyser::modalEnergy ()
// ---------------------------------------------------------------------------
// Print out modal energies per unit area, output by root processor.
// ---------------------------------------------------------------------------
{
  const int_t NC = Geometry::nPert();
  real_t      ek = 0.0, eku = 0.0, ekv = 0.0, ekw = 0.0;
  
  for (int_t i = 0; i < NC; i++) ek += _src -> u[i] -> mode_L2 (0);
 eku=_src -> u[0] -> mode_L2 (0);
 ekv=_src -> u[1] -> mode_L2 (0);
 ekw=_src -> u[2] -> mode_L2 (0);
 
		
   _mdl_strm << setw(10) << _src -> time 
  << setw( 5) << 1 
  << setw(15) << ek
  << setw(15) << eku
  << setw(15) << ekv
  << setw(15) << ekw <<	endl;
}

void Analyser::modalEnergyquad ()
// ---------------------------------------------------------------------------
// Print out modal energies per unit area, output by root processor. Only to the quadrilateral Bachelor vortex.
// ---------------------------------------------------------------------------
{ ofstream outputprofile;
  const int_t NC = Geometry::nPert();
  real_t      ek = 0.0, eku = 0.0, ekv = 0.0, ekw = 0.0;
  int_t i,j;
  //get the mesh array X[x1,x2,...] and Y[y1,y2,...]
  const int_t np    = Geometry::nP();
  const int_t nel    = Geometry::nElmt();
  int_t nX, nY;

  real_t* X;
  real_t* Y;
   X = new real_t [static_cast<size_t>(np*nel)];
   Y = new real_t [static_cast<size_t>(np*nel)];
   _src->u[0]->meshXY( X, Y, nX, nY);
   
   real_t** Ukr;
   Ukr= new real_t* [static_cast<size_t>(NC + 1)];
  for (i=0; i< NC+1; i++)
  Ukr[i] =  new real_t [static_cast<size_t>(nY)];


  
  real_t radius, radiusa=0, radiusb=0, dy;

  real_t dominate_wavenumber=0;
  real_t dominate_energy=0;
  
 //energy in axial mode k
 const int_t N=61;
   real_t axialwaves[N];
   for (i=0;i<N;i++)
   axialwaves[i]=i/20.0;
   real_t energyk[N];
  
   	for(j=0;j<N;j++){
   energyk[j]=0;
	for(i=0;i<NC;i++)
	  energyk[j]+=_src->u[i]->EnergyK(axialwaves[j], X, Y, nX, nY, Ukr[i]);
	if(dominate_energy < energyk[j]) {
	  dominate_energy = energyk[j];  
	  dominate_wavenumber = axialwaves[j];
	  }
	}

// profile of energy with y and radius, Ukr[NC]=(u^2+v^2+w^2)/2/R, her Ukr is different with Ukr in auxfield.c
Veclib::zero (nY,     Ukr[NC], 1);
for(i=0;i<NC;i++){
	_src->u[i]->EnergyK(dominate_wavenumber, X, Y, nX, nY, Ukr[i]);
    Veclib::vvtvp(nY, Ukr[i], 1, Ukr[i], 1, Ukr[NC], 1, Ukr[NC], 1);
}
Veclib::smul (nY, 0.5/(Y[nY-1]-Y[0]) , Ukr[NC], 1, Ukr[NC], 1);

//radius=(integ Ukr*y*dy)/(integ Ukr*dy)
for(i=0;i<nY;i++)
{if (i==0) dy=Y[1]-Y[0];
 else if (i==nY-1) dy=Y[nY-1]-Y[nY-2];
 else dy=(Y[i+1]-Y[i-1])/2;
radiusb+= *(Ukr[NC]+i)*dy;
radiusa+= *(Ukr[NC]+i)*dy*Y[i];
}
radius=radiusa/radiusb;

//output the profile Ukr.
outputprofile.open("energyprofilewithy.txt");
for (i=0;i<nY;i++){
outputprofile << setw(10) << Ukr[NC][i] << setw(10) << Y[i]<<endl;
}
outputprofile.close();
for (i=0; i< NC+1; i++)
delete [] Ukr[i];
delete [] Ukr;


 
  for (int_t i = 0; i < NC; i++)  ek += _src -> u[i] -> mode_L2 (0);
 eku=_src -> u[0] -> mode_L2 (0);
 ekv=_src -> u[1] -> mode_L2 (0);
 ekw=_src -> u[2] -> mode_L2 (0);
 
 
 
  _mdl_strm << setw(10) << _src -> time 
           << setw( 5) << 1 
           << setw(15) << ek<<endl
			    << setw(15) << eku
			    << setw(15) << ekv
					    << setw(15) << ekw <<	endl;
	//for (i=0;i<N;i++)
	//_mdl_strm << setw(15) << energyk[i]
	//	 << endl;
	_mdl_strm << setw(15) << radius << setw(15) << dominate_wavenumber << endl;
	delete [] X;
	delete  [] Y;
	
}


void Analyser::divergence (AuxField** Us) const
// ---------------------------------------------------------------------------
// Print out the velocity field's divergence energy per unit area.  Us
// is used as work area.
// ---------------------------------------------------------------------------
{
  const int_t NC = Geometry::nPert();
  int_t       i;

  if (Geometry::system() == Geometry::Cartesian) {
    for (i = 0; i < NC; i++) {
      *Us[i] = *_src -> u[i];
      Us[i] -> gradient (i);
    }
  } else {
    for (i = 0; i < NC; i++) *Us[i] = *_src -> u[i];
    Us[1] -> mulY();
    for (i = 0; i < NC; i++)  Us[i] -> gradient (i);
    Us[1] -> divY();
    if (NC == 3) Us[2] -> divY();
  }

  if (Geometry::problem() == Geometry::O2_3D_SYMM) *Us[2] *= -1.0;

  for (i = 1; i < NC; i++) *Us[0] += *Us[i];

  cout << "-- Divergence Energy: " << Us[0] -> mode_L2 (0) << endl;
}


void Analyser::estimateCFL () const
// ---------------------------------------------------------------------------
// Estimate and print the peak CFL number, based on zero-mode velocities.
// ---------------------------------------------------------------------------
{
  const real_t CFL_max = 0.7;	// -- Approximate maximum for scheme.
  const real_t SAFETY  = 0.9;	// -- Saftey factor.
  const real_t dt      = Femlib::value ("D_T");
  real_t       CFL_dt, dt_max;
  int_t        percent;

  CFL_dt = max (_src -> u[0] -> CFL (0), _src -> u[1] -> CFL (1));
  if (Geometry::nPert() == 3) CFL_dt = max (CFL_dt, _src -> u[2] -> CFL (2));

  dt_max  = SAFETY * CFL_max / CFL_dt;
  percent = static_cast<int_t>(100.0 * dt / dt_max);

  cout << "-- CFL: "     << CFL_dt * dt;
  cout << ", dt (max): " << dt_max;
  cout << ", dt (set): " << dt;
  cout << " ("           << percent << "%)" << endl;
}
