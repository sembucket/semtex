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

#include <dns.h>

#include "rpo_base.h"
#include "rpo_utils.h"
#include "rpo_preconditioner.h"

#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

static char prog[] = "rpo";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);
void integrate (void (*)(Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
		Domain*, BCmgr*, DNSAnalyser*, FieldForce*);

#define NFIELD 3
#define XMIN 0.0
#define YMIN 0.0
#define YMAX 1.0
#define NSLICE 1

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

void remove_axis(vector<Field*> field, real_t* data_r, real_t* data_i) {
  int plane_j;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = Femlib::ivalue("NELS_X")*elOrd;

  Field::coupleBCs(field[1], field[2], FORWARD);

  for(int field_i = 0; field_i < 3; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
 
      elements_to_logical(Femlib::ivalue("NELS_X"), Femlib::ivalue("NELS_Y"), field[field_i]->plane(plane_i+0), data_r);
      elements_to_logical(Femlib::ivalue("NELS_X"), Femlib::ivalue("NELS_Y"), field[field_i]->plane(plane_i+1), data_i);

      for(int node_x = 0; node_x < nNodesX; node_x++) {
        if(OmitAxis(field_i, 0, plane_j+0)) data_r[node_x] = 0.0;
        if(OmitAxis(field_i, 0, plane_j+1)) data_i[node_x] = 0.0;
      }

      logical_to_elements(Femlib::ivalue("NELS_X"), Femlib::ivalue("NELS_Y"), data_r, field[field_i]->plane(plane_i+0), plane_j, field_i, true);
      logical_to_elements(Femlib::ivalue("NELS_X"), Femlib::ivalue("NELS_Y"), data_i, field[field_i]->plane(plane_i+1), plane_j, field_i, true);
    }
  }

  Field::coupleBCs(field[1], field[2], INVERSE);
}

void rpo_solve(int nSlice, Mesh* mesh, vector<Element*> elmt, BCmgr* bman, Domain* domain, DNSAnalyser* analyst, FieldForce* FF) {
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = Femlib::ivalue("NELS_X")*elOrd;
  int nModesX = nNodesX;
  real_t dx, er, es, ex, ey;
  real_t ckt, skt, cTmp, rTmp;
  const real_t *qx, *wx;
  double dy, dy_sum, norm;
  int_t pt_x, pt_y, el_x, el_y, el_i, el_j, ii, mode_l;
  bool found;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  real_t *data_r, *data_i;
  int np2 = Geometry::nP() * Geometry::nP();
  Vec xi, xj;
  double dz = 2.0*M_PI/360;
  AuxField* uBar;
  vector<AuxField*>   ui;
  vector<AuxField*>   u0;
  vector<AuxField*>   fi;
  int nstep;

  if(!Geometry::procID()) cout << "time step: " << Femlib::value("D_T") << endl;

  context->nField   = NFIELD;
  context->mesh     = mesh;
  context->elmt     = elmt;
  context->domain   = domain;
  context->bman     = bman;
  context->analyst  = analyst;
  context->ff       = FF;
  context->travelling_wave = Femlib::ivalue("TRAV_WAVE");
  context->nElsX    = Femlib::ivalue("NELS_X");
  context->nElsY    = Femlib::ivalue("NELS_Y");

  context->theta_i = new real_t[NSLICE];
  context->phi_i   = new real_t[NSLICE];
  context->tau_i   = new real_t[NSLICE];

  context->nModesX = nNodesX;
  context->xmax    = Femlib::value("XMAX");
  context->beta    = Femlib::ivalue("BETA");

  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    context->theta_i[slice_i] = 0.0;
    context->phi_i[slice_i] = 0.0;
    context->tau_i[slice_i] = Femlib::ivalue("N_STEP") * Femlib::value("D_T");
  }
  context->domain->step++;

  // setup the fourier mapping data
  context->el = new int_t[Femlib::ivalue("NELS_X")*elOrd*Femlib::ivalue("NELS_Y")*elOrd];
  context->r  = new real_t[Femlib::ivalue("NELS_X")*elOrd*Femlib::ivalue("NELS_Y")*elOrd];
  context->s  = new real_t[Femlib::ivalue("NELS_X")*elOrd*Femlib::ivalue("NELS_Y")*elOrd];

  context->rad_weights = new double[context->nElsY*elOrd+1];
  context->rad_coords  = new double[context->nElsY*elOrd+1];

  data_r = new real_t[Femlib::ivalue("NELS_Y")*elOrd*context->nModesX];
  data_i = new real_t[Femlib::ivalue("NELS_Y")*elOrd*context->nModesX];

  dx = (Femlib::value("XMAX")/* - XMIN*/)/nNodesX;
  Femlib::quadrature(&qx, &wx, 0, 0, elOrd+1, GLJ, 0.0, 0.0);

  for(int pt_i = 0; pt_i < Femlib::ivalue("NELS_X")*elOrd*Femlib::ivalue("NELS_Y")*elOrd; pt_i++) {
    pt_x = pt_i%(Femlib::ivalue("NELS_X")*elOrd);
    pt_y = pt_i/(Femlib::ivalue("NELS_X")*elOrd);
    el_x = pt_x/elOrd;
    el_y = pt_y/elOrd;
    el_i = el_y*Femlib::ivalue("NELS_X") + el_x;
    ex = XMIN + pt_x*dx;
    // element y size increases with distance from the boundary
    ey = elmt[el_i]->_ymesh[(pt_y%elOrd)*(elOrd+1)];

    found = false;  
    for(el_j = 0; el_j < mesh->nEl(); el_j++) {
      if(elmt[el_j]->locate(ex, ey, er, es, &work[0], true) && !found) {
        if(fabs(er) < 1.0000000001) {
          context->el[pt_i] = el_j;

          if(er > +0.99999999) er = +0.99999999;
          if(er < -0.99999999) er = -0.99999999;
          if(es > +0.99999999) es = +0.99999999;
          if(es < -0.99999999) es = -0.99999999;

          context->r[pt_i] = er;
          context->s[pt_i] = es;
          found = true;
          //break;
        }
      }
      if(found) break;
    }
    if(!found || el_j==mesh->nEl()) {
      cout << "ERROR! element does not contain point: " << pt_i << "\tx: " << ex << "\ty: " << ey << endl;
      cout << "       pt x: " << pt_x << "\tpt y: " << pt_y << endl;
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

  // add dofs for theta and tau for each time slice
  context->nSlice     = NSLICE;
  context->nField     = NFIELD;
  context->nDofsPlane = nModesX*Femlib::ivalue("NELS_Y")*elOrd;
  context->nDofsSlice = NFIELD * Geometry::nZ() * context->nDofsPlane + 2;
  context->localSize  = nSlice * NFIELD * Geometry::nZProc() * context->nDofsPlane;
  if(!Femlib::ivalue("TRAV_WAVE")) context->nDofsSlice++;
  if(!Geometry::procID()) {
    context->localSize += (nSlice*2);
    if(!Femlib::ivalue("TRAV_WAVE")) context->localSize += nSlice;
  }

  assign_scatter_semtex(context);

/*if(!Geometry::procID()){
int _nr = elOrd*Femlib::ivalue("NELS_Y");
int _nt = Geometry::nZ();
int _nz = elOrd*Femlib::ivalue("NELS_X");
double _dz = Femlib::value("XMAX")/Femlib::ivalue("NELS_X");
char filename[100];
sprintf(filename, "coords.semtex");
ofstream o_file; 
o_file.open(filename);
o_file.precision(16);
for(int iz = 0; iz < _nz; iz++) {
double _z = (iz/elOrd)*_dz + _dz*0.5*(qx[iz%elOrd]+1.0);
for(int it = 0; it < _nt; it++) {
double _theta = ((2.0*M_PI)/Femlib::ivalue("BETA")/Geometry::nZ())*it;
for(int ir = 0; ir < _nr; ir++) {
o_file << context->rad_coords[ir] << "\t" << _theta << "\t" << _z << endl;
}
}
}
o_file.close();
}*/

  // setup the complex fft in the axial direction
  context->data_s = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->data_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->trans_fwd = fftw_plan_dft_1d(context->nModesX, context->data_s, context->data_f, FFTW_FORWARD,  FFTW_ESTIMATE);
  context->trans_bck = fftw_plan_dft_1d(context->nModesX, context->data_f, context->data_s, FFTW_BACKWARD, FFTW_ESTIMATE);

  if(!Geometry::procID()) {
    cout << "nModesX:    " << context->nModesX << endl;
    cout << "nDofsPlane: " << context->nDofsPlane << endl;
  }

  uBar = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), domain->elmt, 'b');
  *uBar = 0.0;
  for(int pl_i = 0; pl_i < Geometry::nZProc(); pl_i++) {
    for(int el_i = 0; el_i < domain->elmt.size(); el_i++) {
      for(int pt_i = 0; pt_i < np2; pt_i++) {
        uBar->plane(pl_i)[el_i*np2+pt_i] = (1.0 - domain->elmt[el_i]->_ymesh[pt_i]*domain->elmt[el_i]->_ymesh[pt_i]);
      }
    }
  }
  uBar->transform(FORWARD);
  context->uBar = uBar;

  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &xi);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &xj);

  ui.resize(3);
  u0.resize(3);
  fi.resize(3);
  for(int field_i = 0; field_i < 3; field_i++) {
    ui[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'U'+field_i);
    fi[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'F'+field_i);
    u0[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'A'+field_i);
    *ui[field_i] = *domain->u[field_i];
    *u0[field_i] = *domain->u[field_i];
  }
  context->ui = ui;
  context->u0 = u0;
  context->fi = fi;

  //base_profile(context, domain->u[0], Femlib::value("BASE_PROFILE_SCALE"), context->uBar);
  *context->ui[0] -= *context->uBar;
  velocity_scales(context);
  *context->ui[0] += *context->uBar;

  cout.precision(12);

  context->c_scale = 1.0;
  RepackX(context, ui, context->theta_i, context->phi_i, context->tau_i, xi, true);
  RPOVecNormL2_Hookstep(context, xi, &norm);
  if(!Geometry::procID()) cout << "|x_0|: " << norm << endl;
  context->c_scale = context->tau_i[0] / norm;

  for(int field_i = 0; field_i < 3; field_i++) *u0[field_i] = *domain->u[field_i];

  // pack the initial data
  RepackX(context, ui, context->theta_i, context->phi_i, context->tau_i, xi, true);
  RPOVecNormL2_Hookstep(context, xi, &norm);
  if(!Geometry::procID()) cout << "initial data norm |x_0|: " << norm << endl;

  // test the final shift
  AuxField::couple(domain->u[1], domain->u[2], INVERSE);
  integrate(skewSymmetric, domain, bman, analyst, FF);
  AuxField::couple(domain->u[1], domain->u[2], FORWARD);
  for(int field_i = 0; field_i < 3; field_i++) *fi[field_i] = *domain->u[field_i];
  phase_shift_x(context, Femlib::value("THETA_0") * (2.0 * M_PI / context->xmax), -1.0, context->domain->u);
  phase_shift_z(context, Femlib::value("PHI_0") * Femlib::ivalue("BETA"), -1.0, context->domain->u);
  for(int field_i = 0; field_i < 3; field_i++) *ui[field_i] -= *domain->u[field_i];
  RepackX(context, ui, context->theta_i, context->phi_i, context->tau_i, xj, true);
  RPOVecNormL2_Hookstep(context, xj, &norm);
  if(!Geometry::procID()) cout << "shifted data norm |x_f - x_0|: " << norm << endl;
  
  // pack the shifted data
/*
  phase_shift_x(context, M_PI, -1.0, domain->u);
  for(int field_i = 0; field_i < 3; field_i++) *ui[field_i] = *domain->u[field_i];
  RepackX(context, ui, context->theta_i, context->phi_i, context->tau_i, xi, true);
  RPOVecNormL2_Hookstep(context, xi, &norm);
  if(!Geometry::procID()) cout << "shifted data norm |x_0|: " << norm << endl;

  if(!Geometry::procID()) cout << "dump to file\n";
  nstep = Femlib::ivalue("N_STEP");
  Femlib::ivalue("N_STEP", 0);
  domain->dump();
  Femlib::ivalue("N_STEP", nstep);
  if(!Geometry::procID()) cout << "done.\n";
*/

  if(Femlib::ivalue("N_STEP") > 1) {
    //AuxField::couple(domain->u[1], domain->u[2], INVERSE);
    //integrate(skewSymmetric, domain, bman, analyst, FF);
    //AuxField::couple(domain->u[1], domain->u[2], FORWARD);
    //for(int field_i = 0; field_i < 3; field_i++) *ui[field_i] = *domain->u[field_i];
    //RepackX(context, ui, context->theta_i, context->phi_i, context->tau_i, xj, true);
    //RPOVecNormL2_Hookstep(context, xj, &norm);
    //cout.precision(12);
    //if(!Geometry::procID()) cout << "theta: " << context->theta_i[0] << "\t|x_f|: " << norm << endl;

    for(int ii = 0; ii <= 360; ii++) {
      for(int field_i = 0; field_i < 3; field_i++) *domain->u[field_i] = *fi[field_i];
      phase_shift_x(context, (1.0*ii)/(2.0*M_PI), -1.0, domain->u);
      for(int field_i = 0; field_i < 3; field_i++) *ui[field_i]  = *u0[field_i];
      for(int field_i = 0; field_i < 3; field_i++) *ui[field_i] -= *domain->u[field_i];
      RepackX(context, ui, context->theta_i, context->phi_i, context->tau_i, xj, true);
      RPOVecNormL2_Hookstep(context, xj, &norm);
      cout.precision(12);
      if(!Geometry::procID()) cout << "theta: " << (1.0*ii)/(2.0*M_PI) << "\t|x_f - x_0|: " << norm << endl;
    }
  }

  delete[] context->el;
  delete[] context->r;
  delete[] context->s;
  delete[] data_r;
  delete[] data_i;
  VecDestroy(&xi);
  VecDestroy(&xj);
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

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);
  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
  domain -> restart ();
  //ROOTONLY domain -> report ();
  
  //Field::coupleBCs(domain->u[1], domain->u[2], FORWARD);
  AuxField::couple(domain->u[1], domain->u[2], FORWARD);

  // solve the newton-rapheson problem
  rpo_solve(NSLICE, mesh, elmt, bman, domain, analyst, FF);
  //delete domain;

  //Field::coupleBCs(domain->u[1], domain->u[2], INVERSE);
  AuxField::couple(domain->u[1], domain->u[2], INVERSE);

  if(!Geometry::procID()) cout << "...done.\n";

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
