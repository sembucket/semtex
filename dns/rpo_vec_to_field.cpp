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
#include "rpo_preconditioner.h"

static char prog[] = "rpo";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);
void integrate (void (*)(Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
		Domain*, BCmgr*, DNSAnalyser*, FieldForce*);

#define TESTING

#define NFIELD 3
#define NSLICE 1
#define THREE 3

static PetscErrorCode RPOVecNormL2_Hookstep(void* ctx,Vec v,PetscScalar* norm) {
  Context* context = (Context*)ctx;
  PetscInt nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  PetscInt ind_i;
  double norm_sq, norm_l_sq = 0.0;
  PetscScalar* vArray;
  Vec vl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl);

  VecScatterBegin(context->global_to_semtex, v, vl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v, vl, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArray(vl, &vArray);
  for(ind_i=0; ind_i<3*nDofsCube_l; ind_i++) {
    norm_l_sq += vArray[ind_i]*vArray[ind_i];
  }
  VecRestoreArray(vl, &vArray);

  norm_sq = 0.0;
  MPI_Allreduce(&norm_l_sq, &norm_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *norm = sqrt(norm_sq);

  VecDestroy(&vl);

  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDot_Hookstep(void* ctx,Vec v1,Vec v2,PetscScalar* dot) {
  Context* context = (Context*)ctx;
  PetscInt nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  PetscInt ind_i;
  double dot_l = 0.0;
  PetscScalar *v1Array, *v2Array;
  Vec vl1, vl2;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl1);
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl2);

  VecScatterBegin(context->global_to_semtex, v1, vl1, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v1, vl1, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterBegin(context->global_to_semtex, v2, vl2, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v2, vl2, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArray(vl1, &v1Array);
  VecGetArray(vl2, &v2Array);
  for(ind_i=0; ind_i<3*nDofsCube_l; ind_i++) {
    dot_l += v1Array[ind_i]*v2Array[ind_i];
  }
  VecRestoreArray(vl1, &v1Array);
  VecRestoreArray(vl2, &v2Array);

  *dot = 0.0;
  MPI_Allreduce(&dot_l, dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  VecDestroy(&vl1);
  VecDestroy(&vl2);

  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDiff_Hookstep(void* ctx,Vec y,Vec F,PetscScalar h) {
  Context* context = (Context*)ctx;
  PetscInt nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  PetscInt ind_i;
  PetscScalar *yArray, *FArray;
  Vec yl, Fl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &yl);
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &Fl);

  VecScatterBegin(context->global_to_semtex, y, yl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, y, yl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterBegin(context->global_to_semtex, F, Fl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, F, Fl, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArray(yl, &yArray);
  VecGetArray(Fl, &FArray);
  for(ind_i=0; ind_i<3*nDofsCube_l; ind_i++) {
    yArray[ind_i] = (yArray[ind_i] - FArray[ind_i])/h;
  }
  VecRestoreArray(yl, &yArray);
  VecRestoreArray(Fl, &FArray);

  VecScatterBegin(context->global_to_semtex, yl, y, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd  (context->global_to_semtex, yl, y, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&yl);
  VecDestroy(&Fl);

  PetscFunctionReturn(0);
}

void rpo_solve(Mesh* mesh, vector<Element*> elmt, BCmgr* bman, FEML* file, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, 
               vector<AuxField*> ui, vector<AuxField*> fi, vector<AuxField*> u0, char* session) {
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = Femlib::ivalue("NELS_X") * elOrd;
  real_t dx, dy, dy_sum, er, es, ex, ey;
  double norm;
  const real_t *qx, *wx;
  int_t pt_x, pt_y, el_x, el_y, el_i, el_j;
  bool found;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  Vec x;
  int np2 = Geometry::nP() * Geometry::nP();
  char filename[100];
  PetscViewer viewer;

  if(!Geometry::procID()) cout << "time step: " << Femlib::value("D_T") << endl;

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
  context->fi       = fi;
  context->u0       = u0;
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
  if( context->travelling_wave && !Geometry::procID()) cout << "      travelling wave solution\n";
  if(!context->travelling_wave && !Geometry::procID()) cout << "NOT a travelling wave solution\n";
  context->iteration = Femlib::ivalue("RPO_LOAD_VEC");

  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    context->theta_i[slice_i] = Femlib::value("THETA_0");
    context->phi_i[slice_i]   = 0.0;
    context->tau_i[slice_i]   = Femlib::ivalue("N_STEP") * Femlib::value("D_T");
  }
  if(!Geometry::procID()) cout << "initial shift: " << context->theta_i[0] * (2.0 * M_PI / context->xmax) << endl;

  if(!context->travelling_wave && !Femlib::ivalue("ITERATIVE")) {
    if(!Geometry::procID()) cout << "ERROR: solver type must be iterative for rpo configuration!!\n";
    abort();
  }

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

  context->uBar = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'p');
  *context->uBar = 0.0;
  for(int pl_i = 0; pl_i < Geometry::nZProc(); pl_i++) {
    for(int el_i = 0; el_i < domain->elmt.size(); el_i++) {
      for(int pt_i = 0; pt_i < np2; pt_i++) {
        context->uBar->plane(pl_i)[el_i*np2+pt_i] = (1.0 - domain->elmt[el_i]->_ymesh[pt_i]*domain->elmt[el_i]->_ymesh[pt_i]);
      }
    }
  }
  context->uBar->transform(FORWARD);

  assign_scatter_semtex(context);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &x);

  // setup the complex fft in the axial direction
  context->data_s = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->data_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->trans_fwd = fftw_plan_dft_1d(context->nModesX, context->data_s, context->data_f, FFTW_FORWARD,  FFTW_ESTIMATE);
  context->trans_bck = fftw_plan_dft_1d(context->nModesX, context->data_f, context->data_s, FFTW_BACKWARD, FFTW_ESTIMATE);

  // dummy initial value
  context->c_scale    = Femlib::value("C_SCALE");
  context->u_scale[0] = Femlib::value("U_SCALE");
  context->u_scale[1] = Femlib::value("V_SCALE");
  context->u_scale[2] = Femlib::value("W_SCALE");

  if(!Geometry::procID()) cout << "loading vectors at iteration: " << context->iteration;

  sprintf(filename, "x_curr_%.4u.vec", context->iteration);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer);
  VecZeroEntries(x);
  VecLoad(x, viewer);
  PetscViewerDestroy(&viewer);

  VecNorm(x, NORM_2, &norm);
  if(!Geometry::procID()) cout << ", |x|: " << norm << endl;

  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x, true);

  VecDestroy(&x);
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
  vector<AuxField*>   fi;  // Solution fields for flow maps at the i time slices
  vector<AuxField*>   u0;  // Initial guess for the solution velocities
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
  fi.resize(NSLICE * THREE);
  u0.resize(NSLICE * THREE);
  delete file;
  delete domain;
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  domain->restart();

  for(int field_i = 0; field_i < THREE/*NFIELD*/; field_i++) {
    ui[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'U'+field_i);
    fi[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'F'+field_i);
    u0[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'A'+field_i);

    *ui[field_i] = *domain->u[field_i];
  }

  // solve the newton-rapheson problem
  rpo_solve(mesh, elmt, bman, file, domain, analyst, FF, ui, fi, u0, session);

  for(int field_i = 0; field_i < THREE; field_i++) *domain->u[field_i] = *ui[field_i];

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
