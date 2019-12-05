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

#include <blas.h>
#include <lapack.h>
#include <mpi.h>
#include <fftw3.h>
#include <dns.h>
#include "rpo_utils.h"

static char prog[] = "rpo";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);
void integrate (void (*)(Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
		Domain*, BCmgr*, DNSAnalyser*, FieldForce*);

#define TESTING 1
#define NFIELD 3
#define THREE 3

#define F77name(x) x ## _

double norm(int nl, double* v) {
  int ii;
  double norm_sq, norm_sq_l = 0.0;

  for(ii = 0; ii < nl; ii++) {
    norm_sq_l += v[ii]*v[ii];
  }
  MPI_Allreduce(&norm_sq_l, &norm_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sqrt(norm_sq);
}

double dot(int nl, double* u, double* v) {
  int ii;
  double dot_l = 0.0, dot_g;

  for(ii = 0; ii < nl; ii++) {
    dot_l += u[ii]*v[ii];
  }
  MPI_Allreduce(&dot_l, &dot_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return dot_g;
}

double normvec(int nl, double* v, double* vn) {
  int ii;
  double norm_q = norm(nl, v);

  for(ii = 0; ii < nl; ii++) {
    vn[ii] = v[ii] / norm_q;
  }
  return norm_q;
}

void mgs(int nk, int nl, double** AK, double** QK) {
  int ii, ki, kj;
  double dot_ij, dot_jj;

  normvec(nl, AK[0], QK[0]);

  for(ki = 1; ki < nk; ki++) {
    for(ii = 0; ii < nl; ii++) QK[ki][ii] = AK[ki][ii];

    for(kj = 0; kj < ki; kj++) {
      dot_ij = dot(nl, QK[ki], QK[kj]);
      dot_jj = dot(nl, QK[kj], QK[kj]); // should be 1!

      for(ii = 0; ii < nl; ii++) QK[ki][ii] -= (dot_ij / dot_jj) * QK[kj][ii];
    }

    normvec(nl, QK[ki], QK[ki]);
  }
}

// Generate an array of orthonormalised vectors, QK, and the Hessenberg, HK
void orthonorm(int nk, int nl, double** AK, double** QK, double** HK) {
  int ii, ki, kj;
  double* v = new double[nl];

  normvec(nl, AK[0], QK[0]);

  for(ki = 0; ki < nk; ki++) {
    for(ii = 0; ii < nl; ii++) v[ii] = AK[ki+1][ii];

    for(kj = 0; kj < ki+1; kj++) {
      HK[kj][ki] = dot(nl, QK[kj], v);

      for(ii = 0; ii < nl; ii++) v[ii] -= HK[kj][ki] * QK[kj][ii];
    }

    HK[ki+1][ki] = normvec(nl, v, QK[ki+1]);
  }

  delete[] v;
}

void repack(Context* context, vector<AuxField*> u, double* v, bool rmv_ubar) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = nex*elOrd;
  int nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int mode_l, index;
  real_t* data_r;
  real_t* data_i;
  AuxField* field;
  double scale;

  data_r = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];
  data_i = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];

  if(rmv_ubar) *u[0] -= *context->uBar;

  for(int field_i = 0; field_i < THREE; field_i++) {
    field = u[field_i];

    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      SEM_to_Fourier(plane_i, context, field, data_r, data_i);

      for(int point_y = 0; point_y < ney*elOrd; point_y++) {
        for(int point_x = 0; point_x < context->nModesX; point_x++) {
          mode_l = (point_x <= context->nModesX/2) ? point_x : point_x - context->nModesX; // fftw ordering of complex data
          scale  = GetScale(context, field_i, plane_i, point_x, point_y);

          index = LocalIndex(context, field_i, plane_i+0, point_x, point_y);
          v[index] = data_r[point_y*context->nModesX+point_x] * scale;
          // divergence free: mean component of radial velocity is 0
          // note that we are in \tilde{} variables, and the nyquist frequency is also 0
          if(field_i  > 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) v[index] = 0.0;

          index = LocalIndex(context, field_i, plane_i+1, point_x, point_y);
          v[index] = data_i[point_y*context->nModesX+point_x] * scale;
          // don't include the nyquist frequency
          if(field_i == 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) v[index] = 0.0;
        }
      }
    }
  }

  if(rmv_ubar) *u[0] += *context->uBar;

  delete[] data_r;
  delete[] data_i;
}

void unpack(Context* context, vector<AuxField*> u, double* v, bool add_ubar) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = nex*elOrd;
  int nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int mode_l, index;
  real_t* data_r;
  real_t* data_i;
  AuxField* field;
  double scale;

  data_r = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];
  data_i = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];

  for(int field_i = 0; field_i < THREE; field_i++) {
    field = u[field_i];

    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      for(int point_y = 0; point_y < ney*elOrd; point_y++) {
        for(int point_x = 0; point_x < context->nModesX; point_x++) {
          mode_l = (point_x <= context->nModesX/2) ? point_x : point_x - context->nModesX; // fftw ordering of complex data
          scale  = GetScale(context, field_i, plane_i, point_x, point_y);

          index = LocalIndex(context, field_i, plane_i+0, point_x, point_y);
          data_r[point_y*context->nModesX+point_x] = v[index] / scale;
          // divergence free: mean component of radial velocity is 0
          // note that we are in \tilde{} variables, and the nyquist frequency is also 0
          if(field_i  > 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) data_r[point_y*context->nModesX+point_x] = 0.0;

          index = LocalIndex(context, field_i, plane_i+1, point_x, point_y);
          data_i[point_y*context->nModesX+point_x] = v[index] / scale;
          // don't include the nyquist frequency
          if(field_i == 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) data_i[point_y*context->nModesX+point_x] = 0.0;
        }
      }

      Fourier_to_SEM(plane_i, context, field, data_r, data_i, field_i);
    }
  }

  if(add_ubar) *u[0] += *context->uBar;

  delete[] data_r;
  delete[] data_i;
}

void residual(Context* context, double* x, double* f) {
  real_t dt;
  int field_i;
  double f_norm, x_norm, dx_norm, runTime, dummy[3], zero = 0.0;
  char filename[100];

  unpack(context, context->ui, x, true);

  // update the starting time for this slice
  context->domain->time = 0.0;
  context->domain->step = 0;
  Femlib::value("t", 0.0);

  // initialise the flow map fields with the solution fields
  for(field_i = 0; field_i < context->nField; field_i++) {
    *context->domain->u[field_i] = *context->ui[field_i];
  }

  AuxField::couple(context->domain->u[1], context->domain->u[2], INVERSE);
#ifdef TESTING
  context->domain->dump();
#endif
  // don't want to call the dns analysis, use custom integrate routine instead
  //delete context->analyst;
  //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
  integrate(convective, context->domain, context->bman, context->analyst, context->ff);
#ifdef TESTING
  context->domain->dump();
#endif
  AuxField::couple(context->domain->u[1], context->domain->u[2], FORWARD);

  phase_shift_x(context, context->theta_i[0] * (2.0 * M_PI / context->xmax), -1.0, context->domain->u);
  phase_shift_z(context, context->phi_i[0], -1.0, context->domain->u);

  // set the residual vector
  for(field_i = 0; field_i < context->nField; field_i++) {
    *context->fi[field_i]  = *context->domain->u[field_i];
    *context->fi[field_i] -= *context->ui[field_i];
  }

  repack(context, context->fi, f, false);

  x_norm = norm(context->localSize, x);
  f_norm = norm(context->localSize, f);
  if(!Geometry::procID()) cout << "\tevaluating function, |x|: " << x_norm << "\t|f|: " << f_norm << endl;
}

void rpo_solve(Mesh* mesh, vector<Element*> elmt, BCmgr* bman, FEML* file, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, 
               vector<AuxField*> ui, vector<AuxField*> fi, vector<AuxField*> u0) {
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = Femlib::ivalue("NELS_X") * elOrd;
  real_t dx, dy, dy_sum, er, es, ex, ey;
  const real_t *qx, *wx;
  int_t pt_x, pt_y, el_x, el_y, el_i, el_j, fd_i, ier, lwork;
  bool found;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  int np2 = Geometry::nP() * Geometry::nP();
  int ki, kj, ji, nk = 20, kdim;
  double **QK, **HK, *AK, *x_o, *x_k, *f_o, *f_k, *wr, *wi, *ZK, *rwork, norm_q, eps = 1.0e-8;

  context->nField   = THREE;
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
  context->nElsX    = Femlib::ivalue("NELS_X");
  context->nElsY    = Femlib::ivalue("NELS_Y");
  context->nModesX = nNodesX;
  context->xmax    = Femlib::value("XMAX");
  if(!Geometry::procID())cout<<"NELS_X: "<<context->nElsX<<", NELS_Y: "<<context->nElsY<<", XMAX: "<<context->xmax<<endl;
  context->theta_i = new real_t[1];
  context->phi_i   = new real_t[1];
  context->tau_i   = new real_t[1];
  context->theta_i[0] = Femlib::value("SHIFT_THETA") / (2.0*M_PI/context->xmax);
  context->phi_i[0]   = Femlib::value("SHIFT_PHI");
  context->tau_i[0]   = Femlib::value("D_T") * Femlib::ivalue("N_STEP");
  if(!Geometry::procID()) cout << "shift theta: " << context->theta_i[0]*(2.0*M_PI/context->xmax)
                               << ", shift phi: " << context->phi_i[0] << ", shift tau: " << context->tau_i[0] << endl;

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
  context->localSize  = THREE * Geometry::nZProc() * context->nDofsPlane;
  MPI_Allreduce(&context->localSize, &context->nDofsSlice, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // setup the complex fft in the axial direction
  context->data_s = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->data_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->trans_fwd = fftw_plan_dft_1d(context->nModesX, context->data_s, context->data_f, FFTW_FORWARD,  FFTW_ESTIMATE);
  context->trans_bck = fftw_plan_dft_1d(context->nModesX, context->data_f, context->data_s, FFTW_BACKWARD, FFTW_ESTIMATE);

  context->uBar = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'b');
  base_profile(context, context->domain->u[0], Femlib::value("BASE_PROFILE_SCALE"), context->uBar);
  *context->ui[0] -= *context->uBar;
  velocity_scales(context);
  *context->ui[0] += *context->uBar;

  /* RUN THE ARNOLDI ITERATION */

  // 0: initialisation
  QK = new double*[nk+1];
  HK = new double*[nk+1];
  for(ki = 0; ki < nk+1; ki++) {
    QK[ki] = new double[context->localSize];
    for(ji = 0; ji < context->localSize; ji++) QK[ki][ji] = 0.0;
    HK[ki] = new double[nk];
    for(kj = 0; kj < nk; kj++) HK[ki][kj] = 0.0;
  }
  x_o = new double[context->localSize];
  x_k = new double[context->localSize];
  f_o = new double[context->localSize];
  f_k = new double[context->localSize];

  // 1: first vector
  repack(context, context->ui, x_o, true);
  residual(context, x_o, f_o);
  normvec(context->localSize, f_o, QK[0]);

  // 2: iterate
  ki = 0;
  do {
    // x + h q_k, h << 1
    for(ji = 0; ji < context->localSize; ji++) x_k[ji] = x_o[ji] + eps * QK[ki][ji];

    // next krylov vector: (f(x + hq) - f(x))/h
    residual(context, x_k, f_k);
    for(ji = 0; ji < context->localSize; ji++) f_k[ji] = (f_k[ji] - f_o[ji]) / eps;

    // orthonormalise (modified gram-schmidt)
    for(ji = 0; ji < context->localSize; ji++) QK[ki+1][ji] = f_k[ji];

    for(kj = 0; kj < ki+1; kj++) {
      HK[kj][ki] = dot(context->localSize, QK[kj], QK[ki+1]);

      for(ji = 0; ji < context->localSize; ji++) QK[ki+1][ji] -= HK[kj][ki] * QK[kj][ji];
    }

    HK[ki+1][ki] = normvec(context->localSize, QK[ki+1], QK[ki+1]);

    if(HK[ki+1][ki] < 1.0e-10) break;
    ki++;
  } while(ki < nk);
  kdim = ki;

  // 3: eigenvalues
  AK = new double[kdim*kdim];
  for(ki = 0; ki < kdim; ki++) {
    for(kj = 0; kj < kdim; kj++) {
      AK[ki*kdim+kj] = HK[ki][kj];
    }
  }
  ZK = new double[kdim*kdim];
  wr = new double[kdim];
  wi = new double[kdim];
  lwork = 4*kdim;
  rwork = new double[lwork];
  
  F77name(dgeev) ("N","V",kdim,AK,kdim,wr,wi,0,1,ZK,kdim,rwork,lwork,ier);
  if(ier && !Geometry::procID()) cout << "Hessenberg eigenvalue error!\n";

  if(!Geometry::procID()) {
    cout << "eigenvalues: \n";
    for(ki = 0; ki < kdim; ki++) {
      cout << "\t" << wr[ki] << "\t+\t" << wi[ki] << "i\n";
    }
  }

  // cleanup
  for(ki = 0; ki < nk+1; ki++) delete[] QK[ki];
  delete[] QK;
  delete[] AK;
  delete[] ZK;
  delete[] x_o;
  delete[] x_k;
  delete[] f_o;
  delete[] f_k;
  delete[] wr;
  delete[] wi;
  delete[] rwork;

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
  char*            fname;
  BoundarySys*     bndry;

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);
  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
  domain -> restart ();
  //ROOTONLY domain -> report ();

  // load in the time slices
  ui.resize(THREE);
  fi.resize(THREE);
  u0.resize(THREE);

  AuxField::couple(domain->u[1], domain->u[2], FORWARD);

  for(int field_i = 0; field_i < THREE; field_i++) {
    ui[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'U'+field_i);
    fi[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'F'+field_i);
    u0[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'A'+field_i);

    *ui[field_i] = *domain->u[field_i];
  }

  rpo_solve(mesh, elmt, bman, file, domain, analyst, FF, ui, fi, u0);
  AuxField::couple(domain->u[1], domain->u[2], INVERSE);

  domain->dump();
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
