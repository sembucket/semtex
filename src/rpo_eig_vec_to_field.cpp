//////////////////////////////////////////////////////////////////////////////
// drive.C: control spectral element DNS for incompressible flows.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// USAGE:
// -----
// dns [options] session
//   options:
//   -h       ... print usage prompt
//   -i       ... use iterative solver for viscous [and pressure] steps
//   -t[t]    ... select time-varying BCs for mode 0 [or all modes]
//   -v[v...] ... increase verbosity level
//   -chk     ... turn off checkpoint field dumps [default: selected]
//   -S|C|N   ... regular skew-symm || convective || Stokes advection
//
// AUTHOR:
// ------
// Hugh Blackburn
// Department of Mechanical & Aerospace Engineering
// Monash University
// Vic 3800
// Australia
// hugh.blackburn@monash.edu
//
// --
// This file is part of Semtex.
//
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

#include <fftw3.h>

#include <dns.h>
#include "rpo_utils.h"

static char prog[] = "rpo";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);
#define NFIELD 3
#define NSLICE 1
#define THREE 3
//#define N_KRYLOV 20
#define N_KRYLOV 27

void rpo_solve(Mesh* mesh, vector<Element*> elmt, BCmgr* bman, FEML* file, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, vector<AuxField*> ui, char* session) {
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = Femlib::ivalue("NELS_X") * elOrd;
  real_t dx, dy, dy_sum, er, es, ex, ey;
  double int_u;
  const real_t *qx, *wx;
  int_t pt_x, pt_y, el_x, el_y, el_i, el_j;
  bool found;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  Vec x, f;
  int np2 = Geometry::nP() * Geometry::nP();
  char filename[100];
  PetscViewer viewer;
  //double EIGVEC[N_KRYLOV] = {-0.37596456, 0.45030068, 0.09739175, 0.01120932, 0.19303401,-0.24779794,
  //0.40186954,-0.29715254, 0.30359398,-0.3578629 , 0.16089847,-0.17532011,
  //0.06274627,-0.10397663, 0.05710285,-0.04256283, 0.02917474,-0.01408484,
 //-0.00115004,-0.00850839};
 double EIGVEC[N_KRYLOV] = {
 0.31125071,-0.05879209, 0.54328523,-0.29586126, 0.0298245,  0.04461417,
  0.09690193,-0.1606904,  0.52707527,-0.2874711,  0.24051219,-0.18644092,
  0.07507548,-0.09542884, 0.06582092,-0.04502774, 0.03914864,-0.00829511,
 -0.0008105 , 0.00330501,-0.00525514, 0.00822228,-0.00696877,-0.00218468,
  0.02817664,-0.0321465,  0.03371018
};

  context->nSlice   = NSLICE;
  context->nField   = NFIELD;
  context->mesh     = mesh;
  context->elmt     = elmt;
  context->domain   = domain;
  context->bman     = bman;
  context->file     = file;
  context->analyst  = analyst;
  context->ff       = FF;
  context->ui       = ui;
  context->travelling_wave = Femlib::ivalue("TRAV_WAVE");
  context->build_dx = false;
  context->nElsX    = Femlib::ivalue("NELS_X");
  context->nElsY    = Femlib::ivalue("NELS_Y");
  context->beta     = Femlib::ivalue("BETA");

  context->theta_i = new real_t[NSLICE];
  context->phi_i   = new real_t[NSLICE];
  context->tau_i   = new real_t[NSLICE];

  context->nModesX = nNodesX;
  context->xmax    = Femlib::value("XMAX");
  if(!Geometry::procID())cout<<"NELS_X: "<<context->nElsX<<", NELS_Y: "<<context->nElsY<<", XMAX: "<<context->xmax<<endl;

  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    context->theta_i[slice_i] = Femlib::value("THETA_0");
    context->phi_i[slice_i]   = 0.0;
    context->tau_i[slice_i]   = Femlib::ivalue("N_STEP") * Femlib::value("D_T");
  }
  if(!Geometry::procID()) cout << "initial shift: " << context->theta_i[0] * (2.0 * M_PI / context->xmax) << endl;

  // setup the fourier mapping data
  context->el = new int_t[context->nElsX*elOrd*context->nElsY*elOrd];
  context->r  = new real_t[context->nElsX*elOrd*context->nElsY*elOrd];
  context->s  = new real_t[context->nElsX*elOrd*context->nElsY*elOrd];

  context->rad_weights = new double[context->nElsY*elOrd+1];
  context->rad_coords  = new double[context->nElsY*elOrd+1];

  dx = (context->xmax/* - XMIN*/)/nNodesX;
  Femlib::quadrature(&qx, &wx, 0, 0, elOrd+1, GLJ, 0.0, 0.0);

  for(int pt_i = 0; pt_i < context->nElsX*elOrd*context->nElsY*elOrd; pt_i++) {
    pt_x = pt_i%(context->nElsX*elOrd);
    pt_y = pt_i/(context->nElsX*elOrd);
    el_x = pt_x/elOrd;
    el_y = pt_y/elOrd;
    el_i = el_y*context->nElsX + el_x;
    ex = /*XMIN +*/ pt_x*dx;
    // element y size increases with distance from the boundary
    ey = elmt[el_i]->_ymesh[(pt_y%elOrd)*(elOrd+1)];

    //ex += 1.0e-8; ey += 1.0e-8;
    found = false;  
    for(el_j = 0; el_j < mesh->nEl(); el_j++) {
      er = es = 0.0;
      if(elmt[el_j]->locate(ex, ey, er, es, &work[0], true) && !found) {
        if(fabs(er) < 1.0000000001 && fabs(es) < 1.0000000001) {
          context->el[pt_i] = el_j;

          if(er > +0.99999999) er = +0.99999999;
          if(er < -0.99999999) er = -0.99999999;
          if(es > +0.99999999) es = +0.99999999;
          if(es < -0.99999999) es = -0.99999999;

          context->r[pt_i] = er;
          context->s[pt_i] = es;
          found = true;
        }
      }
      if(found) break;
    }
    if(!found && !Geometry::procID()) {
      cout << Geometry::procID() << "\t:ERROR! element does not contain point: " << pt_i << "\tx: " << ex << "\ty: " << ey << endl;
      cout << "\tfound = " << found << endl;
      cout << "\tel i:   " << el_j << "\tnum els: " << mesh->nEl() << endl;
      cout << "\tpt x: " << pt_x << "\tpt y: " << pt_y << endl;
      abort();
    }
  }

  // compute the radial weights
  for(int pt_i = 0; pt_i <= context->nElsY*elOrd; pt_i++) context->rad_weights[pt_i] = 0.0;

  dy_sum = 0.0;
  for(int el_y = 0; el_y < context->nElsY; el_y++) {
    el_i = el_y*context->nElsX;
    dy = fabs(elmt[el_i]->_ymesh[elOrd*(elOrd+1)] - elmt[el_i]->_ymesh[0]);
    if(!Geometry::procID())printf("%d\tdy: %g\n",el_y,dy);
    for(int qp_i = 0; qp_i <= elOrd; qp_i++) {
      pt_y = el_y*elOrd + qp_i;
      context->rad_weights[pt_y] += 0.5 * dy * wx[qp_i];
      context->rad_coords[pt_y]   = dy_sum + 0.5 * dy * (qx[qp_i]+1.0);
    }
    dy_sum += dy;
  }
  if(!Geometry::procID())printf("\tdy_sum: %g\n",dy_sum);

  // add dofs for theta and tau for each time slice
  context->nDofsPlane = context->nModesX*context->nElsY*elOrd;
  context->nDofsSlice = context->nField * Geometry::nZ() * context->nDofsPlane + 1; // add azimuthal phase shift dof
  context->nDofsSlice++; // add axial phase shift dof

  if(!context->travelling_wave) context->nDofsSlice++; // add temporal phase shift dof

  context->localSize  = NSLICE * NFIELD * Geometry::nZProc() * context->nDofsPlane;

  if(!Geometry::procID()) context->localSize += NSLICE; // azimuthal phase shift dof
  if(!Geometry::procID()) context->localSize += NSLICE; // axial phase shift dof
  if(!context->travelling_wave && !Geometry::procID()) context->localSize += NSLICE; // temporal phase shift dof

  assign_scatter_semtex(context);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &x);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &f);

  // setup the complex fft in the axial direction
  context->data_s = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->data_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->trans_fwd = fftw_plan_dft_1d(context->nModesX, context->data_s, context->data_f, FFTW_FORWARD,  FFTW_ESTIMATE);
  context->trans_bck = fftw_plan_dft_1d(context->nModesX, context->data_f, context->data_s, FFTW_BACKWARD, FFTW_ESTIMATE);

  // velocity component scales
  context->c_scale    = Femlib::value("C_SCALE");
  context->u_scale[0] = Femlib::value("U_SCALE");
  context->u_scale[1] = Femlib::value("V_SCALE");
  context->u_scale[2] = Femlib::value("W_SCALE");

  // load state from petsc vectors
  if(!Geometry::procID()) cout << "loading krylov vectors: " << N_KRYLOV << endl;

  VecZeroEntries(f);
  for(int vec_i = 0; vec_i < N_KRYLOV; vec_i++) {
    sprintf(filename, "krylov_%.4u.vec", vec_i);
    if(!Geometry::procID()) cout << "loading krylov vector: " << filename << ".....\t";
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer);
    VecZeroEntries(x);
    VecLoad(x, viewer);
    VecNorm(x, NORM_2, &int_u);
    if(!Geometry::procID()) cout << "val: " << EIGVEC[vec_i] << ",\t|x|: " << int_u;
    VecAXPY(f, EIGVEC[vec_i], x);
    VecNorm(f, NORM_2, &int_u);
    if(!Geometry::procID()) cout << ",\t|f|: " << int_u << endl;
    PetscViewerDestroy(&viewer);
  }
  RepackX(context, context->ui, &context->f_theta, &context->f_phi, &context->f_tau, f, false);

  for(int fld_i = 0; fld_i < 3; fld_i++) {
    context->ui[fld_i]->transform(INVERSE);

    *context->domain->u[3] = 0.0;
    context->domain->u[3]->times(*context->ui[fld_i], *context->ui[fld_i]);
    context->domain->u[3]->transform(FORWARD);
    int_u = context->domain->u[3]->integral();
    if(!Geometry::procID()) cout << "int u " << fld_i << ":\t" << int_u << endl;

    context->ui[fld_i]->transform(FORWARD);
    context->ui[fld_i]->zeroNyquist();
  }

  VecDestroy(&x);
  VecDestroy(&f);
  delete[] context->el;
  delete[] context->r;
  delete[] context->s;
  delete[] context->rad_weights;
  delete[] context->rad_coords;
}

int main (int argc, char** argv) {
#ifdef _GNU_SOURCE
  feenableexcept (FE_OVERFLOW);    // -- Force SIG8 crash on FP overflow.
#endif
  char*            session;
  bool             freeze = false;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  Domain*          domain;
  DNSAnalyser*     analyst;
  FieldForce*      FF;
  static char      help[] = "petsc";
  vector<AuxField*>   ui;  // Solution fields for velocities, pressure at the i time slices
  char             session_i[100];
  FEML*            file_i;
  char*            fname;
  BoundarySys*     bndry;

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);
  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
  //domain -> restart ();
  //ROOTONLY domain -> report ();

  // load in the time slices
  ui.resize(NSLICE * THREE);
  delete file;
  delete domain;
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  domain->restart();
  AuxField::couple(domain->u[1], domain->u[2], FORWARD);

  for(int field_i = 0; field_i < THREE/*NFIELD*/; field_i++) {
    ui[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'U'+field_i);
  }

  // solve the newton-rapheson problem
  rpo_solve(mesh, elmt, bman, file, domain, analyst, FF, ui, session);

  for(int field_i = 0; field_i < THREE; field_i++) {
    *domain->u[field_i] = *ui[field_i];
  }
  AuxField::couple(domain->u[1], domain->u[2], INVERSE);
  domain->dump();
  delete file_i;
  delete domain;

  PetscFinalize();

  return EXIT_SUCCESS;
}

static void getargs (int    argc   ,
		     char** argv   ,
		     bool&  freeze ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char       buf[StrMax];
  const char routine[] = "getargs";
  const char usage[]   = "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -f       ... freeze velocity field (scalar advection/diffusion only)\n"
    "  -i       ... use iterative solver for viscous steps\n"
    "  -v[v...] ... increase verbosity level\n"
    "  -chk     ... turn off checkpoint field dumps [default: selected]\n"
    "  -S|C|N   ... regular skew-symm || convective || Stokes advection\n";

  Femlib::ivalue ("ADVECTION", 1); // -- Default is alternating skew symmetric.

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'S': Femlib::ivalue ("ADVECTION", 0); break;
    case 'C': Femlib::ivalue ("ADVECTION", 2); break;
    case 'N': Femlib::ivalue ("ADVECTION", 3); break;
    case 'f':
      freeze = true;
      break;
    case 'i':
      do			// -- Only allowing ITERATIVE=1 (Viscous).
	Femlib::ivalue ("ITERATIVE", 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE",   Femlib::ivalue ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv))     Femlib::ivalue ("CHKPOINT",    0);
      else { fprintf (stdout, usage, prog); exit (EXIT_FAILURE); }
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  //if   (argc != 1) message (routine, "no session definition file", ERROR);
  //else             session = *argv;
  session = *argv;

  Femlib::value ("DTBDX", 0.0);
}

static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			Domain*&          domain ,
			FieldForce*&      FF     )
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed above in order of creation.
// ---------------------------------------------------------------------------
{

  const char routine[] = "preprocess";
  const int_t        verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  int_t              i, np, nz, nel, procid, seed;

  // -- Initialise problem and set up mesh geometry.
  VERBOSE cout << "Building mesh ..." << endl;
  file = new FEML (session);
  mesh = new Mesh (file);
  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.
  VERBOSE cout << "Setting geometry ... ";
  nel   =  mesh -> nEl();
  np    =  Femlib::ivalue ("N_P");
  nz    =  Femlib::ivalue ("N_Z");
  space = (Femlib::ivalue ("CYLINDRICAL")) ?
    Geometry::Cylindrical : Geometry::Cartesian;

  Geometry::set (np, nz, nel, space);
  VERBOSE cout << "done" << endl;

  // -- If token RANSEED > 0 then initialize the random number
  //    generator based on wall clock time and process ID (i.e. a "truly"
  //    pseudo-random number).  NB: it is important to have done this
  //    before any other possible call to random number routines.

  if (Femlib::ivalue("RANSEED") > 0) {
    procid = Geometry::procID();
    seed   = -abs((procid + 1) * (char) time(NULL));
  } else seed = -1;
  Veclib::ranInit (seed);

  // -- Build all the elements.
  VERBOSE cout << "Building elements ... ";
  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);
  VERBOSE cout << "done" << endl;

  // -- Build all the boundary condition applicators.
  VERBOSE cout << "Building boundary condition manager ..." << endl;
  bman = new BCmgr (file, elmt);
  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.
  VERBOSE cout << "Building domain ..." << endl;
  domain = new Domain (file, elmt, bman);
  VERBOSE cout << "done" << endl;

  // -- Build field force.
  VERBOSE cout << "Building field force ..." << endl;
  FF = new FieldForce (domain, file);
  VERBOSE cout << "done" << endl;
}
