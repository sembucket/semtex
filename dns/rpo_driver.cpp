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

#define X_FOURIER
#define NFIELD 3
#define XMIN 0.0
#define XMAX (2.0*M_PI)
#define YMIN 0.0
#define YMAX 1.0
#define NELS_X 30
#define NELS_Y 7
//#define NSLICE 16
#define NSLICE 1
//#define NSTEPS 3200
#define NSTEPS ((40*16)/NSLICE)

void build_constraints(Context* context, Vec x_delta, double* f_theta, double* f_phi, double* f_tau) {
  int          elOrd       = Geometry::nP() - 1;
  int          nNodesX     = NELS_X*elOrd;
  int          nModesX     = nNodesX/2;// + 2;
  int          slice_i     = 0;
  int          index;
  int          nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int          nl          = context->nField * nDofsCube_l;
  int          plane_j;
  double       k_x, k_z;
  double       f_theta_l, f_phi_l, f_tau_l;
  double*      rx          = new double[nl];
  double*      rz          = new double[nl];
  double*      rt          = new double[nl];
  double*      data_r      = new double[context->nDofsPlane];
  double*      data_i      = new double[context->nDofsPlane];
  PetscScalar* xArray;
  Vec          xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecGetArray(xl, &xArray);

  for(int field_i = 0; field_i < context->nField; field_i++) {
    *context->domain->u[field_i] = *context->ui[slice_i * context->nField + field_i];
  }
  Femlib::ivalue("N_STEP", 1);

  delete context->analyst;
  context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
  integrate(skewSymmetric, context->domain, context->bman, context->analyst, context->ff);
  for(int field_i = 0; field_i < context->nField; field_i++) {
    *context->domain->u[field_i] -= *context->ui[slice_i * context->nField + field_i];
    *context->domain->u[field_i] *= 1.0 / Femlib::value("D_T");
  } 

  for(int dof_i = 0; dof_i < nl; dof_i++) { rx[dof_i] = rz[dof_i] = rt[dof_i] = 0.0; }

  for(int field_i = 0; field_i < context->nField; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      SEM_to_Fourier(plane_i, context, context->ui[slice_i*context->nField+field_i], data_r, nModesX);
      for(int node_j = 0; node_j < NELS_Y * elOrd; node_j++) {
        for(int mode_i = 0; mode_i < nModesX; mode_i++) {
          index = field_i * nDofsCube_l + plane_i * context->nDofsPlane + node_j * nModesX + mode_i;
          k_x = context->theta_i[slice_i] * (mode_i / 2);
          if(mode_i % 2 == 0) {
            rx[index] = -k_x * data_r[node_j * nModesX + mode_i + 1];
          } else {
            rx[index] = +k_x * data_r[node_j * nModesX + mode_i - 1];
          }
        }
      }
    }
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i+=2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
      k_z = context->phi_i[slice_i] * (plane_j / 2);
      SEM_to_Fourier(plane_i+0, context, context->ui[slice_i*context->nField+field_i], data_r, nModesX);
      SEM_to_Fourier(plane_i+1, context, context->ui[slice_i*context->nField+field_i], data_i, nModesX);
      for(int dof_i = 0; dof_i < NELS_Y * elOrd * nModesX; dof_i++) {
        index = field_i * nDofsCube_l + (plane_i+0) * context->nDofsPlane + dof_i;
        rz[index] = -k_z * data_i[dof_i];
        index = field_i * nDofsCube_l + (plane_i+1) * context->nDofsPlane + dof_i;
        rz[index] = +k_z * data_r[dof_i];
      }
    }
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      SEM_to_Fourier(plane_i, context, context->domain->u[field_i], data_r, nModesX);
      for(int dof_i = 0; dof_i < context->nDofsPlane; dof_i++) {
        index = field_i * nDofsCube_l + plane_i * context->nDofsPlane + dof_i;
        rt[index] = data_r[dof_i];
      }
    }
  }

  f_theta_l = f_phi_l = f_tau_l = 0.0;
  for(int dof_i = 0; dof_i < nl; dof_i++) {
    f_theta_l += rx[dof_i] * xArray[dof_i];
    f_phi_l   += rz[dof_i] * xArray[dof_i];
    f_tau_l   += rt[dof_i] * xArray[dof_i];
  }
  MPI_Allreduce(&f_theta_l, f_theta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&f_phi_l,   f_phi,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&f_tau_l,   f_tau,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  VecRestoreArray(xl, &xArray);

  delete[] rx;
  delete[] rz;
  delete[] rt;
  delete[] data_r;
  delete[] data_i;
  VecDestroy(&xl);
}

PetscErrorCode _snes_jacobian(SNES snes, Vec x, Mat J, Mat P, void* ctx) {
  Context* context = (Context*)ctx;

  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  J, MAT_FINAL_ASSEMBLY);

  if(context->build_PC) {
    build_preconditioner_ffs(context, P);
    //context->build_PC = false;
  }
  return 0;
}

PetscErrorCode _snes_function(SNES snes, Vec x, Vec f, void* ctx) {
  Context* context = (Context*)ctx;
  real_t dt;
  int slice_i, slice_j, field_i, mode_i, mode_j, dof_i, nStep;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2;// + 2;
  real_t time = 0.0;
  real_t f_norm, x_norm;
  double f_theta[1], f_phi[1], f_tau[1];
  Vec x_delta;
  PetscScalar *xArray;

  // create the \delta x vector for use in the assembly of the constraints into the residual vector
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nSlice * context->nDofsSlice, &x_delta);
  VecCopy(x, x_delta);
  VecAXPY(x_delta, -1.0, context->x_prev);
  {
    double delta_x_norm;
    VecNorm(x_delta, NORM_2, &delta_x_norm);
    if(!Geometry::procID()) cout << "|delta x|: " << delta_x_norm << endl;
  }
  VecCopy(x, context->x_prev);

  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);

  Femlib::value("D_T", 0.02); // 80x simulation value
  dt = Femlib::value ("D_T");

  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    // update the starting time for this slice
    context->domain->time = time;
    Femlib::value("t", time);

    // initialise the flow map fields with the solution fields
    for(field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] = *context->ui[slice_i * context->nField + field_i];
    }

    // solve the flow map for time tau_i
    MPI_Bcast(&context->tau_i[slice_i]  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&context->theta_i[slice_i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&context->phi_i[slice_i]  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    nStep = (int)(context->tau_i[slice_i]/dt);
    Femlib::ivalue("N_STEP", nStep);
    Femlib::value("D_T", context->tau_i[slice_i]/nStep);
    context->domain->step = 0;

    if(!Geometry::procID()) {
      cout << "time: " << time << "\tslice: " << slice_i << "\tnstep: " << Femlib::ivalue("N_STEP") << endl;
      cout << scientific << "\ttau:   " << context->tau_i[slice_i]   
                         << "\ttheta: " << context->theta_i[slice_i] 
                         << "\tphi:   " << context->phi_i[slice_i] << endl;
    }

    // don't want to call the dns analysis, use custom integrate routine instead
    delete context->analyst;
    context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
    integrate(skewSymmetric, context->domain, context->bman, context->analyst, context->ff);

#ifdef X_FOURIER
    // phase shift in theta (axial direction)
    phase_shift_x(context, context->theta_i[slice_i], -1.0, context->domain->u);
#endif
    // phase shift in phi (azimuthal direction)
    phase_shift_z(context, context->phi_i[slice_i],   -1.0, context->domain->u);

    // set f
    for(field_i = 0; field_i < context->nField; field_i++) {
      slice_j = (slice_i + 1)%context->nSlice;
      *context->fi[slice_i * context->nField + field_i]  = *context->domain->u[field_i];
      *context->fi[slice_i * context->nField + field_i] -= *context->ui[slice_j * context->nField + field_i];
    }

    time += context->tau_i[slice_i];
  }

  build_constraints(context, x_delta, f_theta, f_phi, f_tau);

  RepackX(context, context->fi, f_theta, f_phi, f_tau, f);
  Femlib::value("D_T", dt);

#ifdef TESTING
  char session_i[100];
  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    sprintf(session_i, "tube8_tw_rpo_%u", slice_i + 1);
    FEML* file_i = new FEML(session_i);
    Domain* dom = new Domain(file_i, context->elmt, context->bman);
    for(field_i = 0; field_i < context->nField; field_i++) {
      dom->u[field_i] = context->ui[slice_i*context->nField+field_i];
    }
    dom->dump();
    delete file_i;
    delete dom;
  }
#endif
 
  VecNorm(x, NORM_2, &x_norm);
  VecNorm(f, NORM_2, &f_norm);
  if(!Geometry::procID()) cout << "evaluating function, |x|: " << x_norm << "\t|f|: " << f_norm << endl;

  VecDestroy(&x_delta);

  return 0;
}

void rpo_solve(int nSlice, Mesh* mesh, vector<Element*> elmt, BCmgr* bman, FEML* file, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, 
               vector<Field*> ui, vector<Field*> fi, vector<Field*> uj)
{
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2;// + 2;
  real_t dx, er, es, ex, ey;
  const real_t* qx;
  int_t pt_x, pt_y, el_x, el_y, el_i, el_j;
  const bool guess = false;
  bool found;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  Vec x, f;
  Mat P;
  SNES snes;
  IS *is_s, *is_u, *is_p;

  if(!Geometry::procID()) cout << "time step: " << Femlib::value("D_T") << endl;

  context->nSlice   = nSlice;
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
  context->uj       = uj;
  context->build_PC = true;
  context->is_s     = NULL;
  context->is_u     = NULL;
  context->is_p     = NULL;

  context->theta_i = new real_t[nSlice];
  context->phi_i   = new real_t[nSlice];
  context->tau_i   = new real_t[nSlice];

  Femlib::value("D_T", 0.02); // 80x simulation value
  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
    context->theta_i[slice_i] = 0.0;
    context->phi_i[slice_i] = 0.0;
    context->tau_i[slice_i] = NSTEPS*Femlib::value("D_T");
    if(!Geometry::procID()) cout << "slice: " << slice_i << "\ttau: " << context->tau_i[slice_i] << endl;
  }

  // setup the fourier mapping data
  context->el = new int_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->r  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->s  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];

  dx = (XMAX - XMIN)/nNodesX;
  Femlib::quadrature(&qx, 0, 0, 0, elOrd+1, GLJ, 0.0, 0.0);

  for(int pt_i = 0; pt_i < NELS_X*elOrd*NELS_Y*elOrd; pt_i++) {
    pt_x = pt_i%(NELS_X*elOrd);
    pt_y = pt_i/(NELS_X*elOrd);
    el_x = pt_x/elOrd;
    el_y = pt_y/elOrd;
    el_i = el_y*NELS_X + el_x;
    ex = XMIN + pt_x*dx;
    // element y size increases with distance from the boundary
    ey = elmt[el_i]->_ymesh[(pt_y%elOrd)*(elOrd+1)];

    found = false;  
    for(el_j = 0; el_j < mesh->nEl(); el_j++) {
      // pass er and es by reference?
      if(elmt[el_j]->locate(ex, ey, er, es, &work[0], guess) && !found) {
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

  // add dofs for theta and tau for each time slice
#ifdef X_FOURIER
  context->nDofsPlane = nModesX*NELS_Y*elOrd;
  context->nDofsSlice = context->nField * Geometry::nZ() * context->nDofsPlane + 3;
  context->localSize  = nSlice * NFIELD * Geometry::nZProc() * context->nDofsPlane;
  if(!Geometry::procID()) context->localSize += (nSlice * 3);
#else
  context->nDofsPlane = NELS_X*elOrd*NELS_Y*elOrd;
  context->nDofsSlice = context->nField * Geometry::nZ() * context->nDofsPlane + 2;
  context->localSize  = nSlice * NFIELD * Geometry::nZProc() * context->nDofsPlane;
  if(!Geometry::procID()) context->localSize += (nSlice * 2);
#endif

  assign_scatter_semtex(context);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, nSlice * context->nDofsSlice, &x);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, nSlice * context->nDofsSlice, &f);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, nSlice * context->nDofsSlice, &context->x_prev);

  MatCreate(MPI_COMM_WORLD, &P);
  MatSetType(P, MATMPIAIJ);
  MatSetSizes(P, context->localSize, context->localSize, nSlice * context->nDofsSlice, nSlice * context->nDofsSlice);
  MatMPIAIJSetPreallocation(P, 2*nModesX, PETSC_NULL, 2*nModesX, PETSC_NULL);

  SNESCreate(MPI_COMM_WORLD, &snes);
  SNESSetFunction(snes, f,    _snes_function, (void*)context);
  SNESSetJacobian(snes, P, P, _snes_jacobian, (void*)context);
  SNESSetType(snes, SNESNEWTONTR);
  SNESSetNPCSide(snes, PC_LEFT);
  SNESSetFromOptions(snes);

  context->snes = snes;
  context->x_norm = 100.0;
  RepackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);
  //VecNorm(x, NORM_2, &context->x_norm);
  //if(!Geometry::procID()) cout << "|x_0|: " << context->x_norm << endl;
  VecCopy(x, context->x_prev);
  SNESSolve(snes, NULL, x);
  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);

  VecDestroy(&x);
  VecDestroy(&f);
  MatDestroy(&P);
  VecDestroy(&context->x_prev);
  delete[] context->el;
  delete[] context->r;
  delete[] context->s;
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
  vector<Field*>   ui;  // Solution fields for velocities, pressure at the i time slices
  vector<Field*>   fi;  // Solution fields for flow maps at the i time slices
  vector<Field*>   uj;  // Dummy fields for updating the solution fields within the rpo solver
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
  ui.resize(NSLICE * NFIELD);
  fi.resize(NSLICE * NFIELD);
  uj.resize(NSLICE * NFIELD);
  delete file;
  delete domain;
  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    //sprintf(session_i, "%s.%u", session, slice_i + 1);
    sprintf(session_i, "%s.%u", session, 4*slice_i);
    file_i = new FEML(session_i);
    domain = new Domain(file_i, elmt, bman);
    domain->restart();
    for(int field_i = 0; field_i < NFIELD; field_i++) {
      ui[slice_i*NFIELD+field_i] = domain->u[field_i];

      strcpy ((fname = new char [strlen (bman -> field()) + 1]), bman -> field());

      bndry = new BoundarySys(bman, elmt, fname[0]);
      fi[slice_i*NFIELD+field_i] = new Field(bndry, new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, fname[0]);

      bndry = new BoundarySys(bman, elmt, fname[0]);
      uj[slice_i*NFIELD+field_i] = new Field(bndry, new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, fname[0]);
    }
    delete file_i;
    delete domain;
  }

  // allocate the temporary fields
  sprintf(session_i, "%s.0", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  domain->restart();

  // solve the newton-rapheson problem
  rpo_solve(NSLICE, mesh, elmt, bman, file, domain, analyst, FF, ui, fi, uj);
  delete file_i;
  delete domain;

  // dump the output
  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    sprintf(session_i, "%s_rpo_%u", session, 4*slice_i);
    file_i = new FEML(session_i);
    domain = new Domain(file_i, elmt, bman);
    for(int field_i = 0; field_i < NFIELD; field_i++) {
      domain->u[field_i] = ui[slice_i*NFIELD+field_i];
    }
    domain->dump();
    delete file_i;
    delete domain;
  } 

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
