///////////////////////////////////////////////////////////////////////////////
// integrate.cpp: Unsteady Navier--Stokes solver, using
// "stiffly-stable" time integration [1,2].  Geometries may be 2- or
// 3-dimensional, Cartesian or cylindrical [3].  Fourier expansions
// are used in the homogeneous (z) direction.  This file provides
// integrate as a call-back routine; after initialisation, integrate
// may be called repeatedly without reinitialising internal storage.
//
// For cylindrical coordinates (Fourier in azimuth):
//   u <==> axial     velocity,  x <==> axial     coordinate direction,
//   v <==> radial    velocity,  y <==> radial    coordinate direction,
//   w <==> azimuthal velocity,  z <==> azimuthal coordinate direction.
//
// For Cartesian coordinates (Fourier in z):
//   u <==> x-component  velocity
//   v <==> y-component  velocity
//   w <==> z-component  velocity
//
// In either system, the w velocity component is optional for 2D
// (N_Z=1) (i.e. can have 2D2C or 2D3C).  If 3D (N_Z > 1), w should
// appear in session.
//
// Optionally integrate concentration of advected scalar field c.
//
// Copyright (c) 1994 <--> $Date: 2019/06/21 13:22:31 $, Hugh Blackburn
//
// REFERENCES
// ----------
// [1] Karniadakis, Israeli & Orszag (1991) "High-order splitting methods
//     for the incompressible Navier--Stokes equations", JCP 97:414--443
// [2] Guermond & Shen (2003) "Velocity correction projection methods for
//     incompressible flows", SIAM J Numer Anal 41:112-134
// [3] Blackburn & Sherwin (2004) "Formulation of a Galerkin spectral
//     element--Fourier method for three-dimensional incompressible flows
//     in cylindrical geometries", JCP 179:759-778
// [4] Dong, Karniakadis & Chryssostomides (2014) "A robust and
//     accurate outflow boundary condition for incompressible flow
//     simulations on severely-truncated unbounded domains", JCP 261:83-105.
// [5] Blackburn, Lee, Albrecht & Singh (2019) "Semtex: a spectral
//     element–Fourier solver for the incompressible Navier–Stokes
//     equations in cylindrical or Cartesian coordinates", CPC.
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
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: integrate.cpp,v 9.2 2019/06/21 13:22:31 hmb Exp $";

#include <dns.h>
#include <mpi.h>
#include "rpo_base.h"

typedef ModalMatrixSys Msys;

// -- File-scope constants and routines:

static int_t NDIM, NCOM, NORD, NADV, PIND, CIND;
static bool  C3D;

static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const int_t, AuxField*, Msys*);

bool alloc_rk = true;
bool alloc_pc = true;
AuxField**        du;    // 9 component tensor for the velocity gradients
AuxField**        vort;  // 3 component vector for the vorticity
AuxField**        visc;  // 3 component vector for the explicit viscous terms
AuxField**        rhs;   // 3 component vector for the rhs
AuxField**        phi;   // multi-stage pressure storage
AuxField**        phi_c; // multi-stage pressure storage (with correction)
AuxField**        v_i;   // multi-stage predictor substep storage
AuxField**        u_i;   // multi-stage corrector substep storage
AuxField**        F_i;   // multi-stage forcing term storage
AuxField**        u_o;   // initial velocities
vector<AuxField*> __int_ke_vec;
AuxField*         __int_ke_tmp;
AuxField**        __int_err_prev;
AuxField*         __int_err_temp;
AuxField*         __divg_rhs_tmp;
AuxField**        __visc_rhs_tmp;
AuxField**        __diag_vel_tmp;
AuxField*         tmp;   // temporary field for intermediate evaluations
AuxField**        __int_nf_tmp = NULL;
vector<AuxField*> vel0;

void init_fields(Domain* domain) {
  const int_t   nmodes = Geometry::nModeProc();
  const int_t   base   = Geometry::baseMode();
  const real_t  beta   = Femlib::value("BETA");
  const int_t   ntime  = Femlib::ivalue("N_TIME");

  if(!alloc_rk) return;

  if(!Geometry::procID()) cout << "allocating rk integration fields....\n";

  du  = new AuxField*[9];
  F_i = new AuxField*[9];
  u_i = new AuxField*[9];
  v_i = new AuxField*[9];
  for(int ii = 0; ii < 9; ii++) {
    du[ii]  = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    F_i[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    u_i[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    v_i[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
  }
  vort  = new AuxField*[3];
  visc  = new AuxField*[3];
  rhs   = new AuxField*[3];
  phi   = new AuxField*[3];
  phi_c = new AuxField*[3];
  u_o   = new AuxField*[3];
  __int_err_prev = new AuxField*[3];
  __visc_rhs_tmp = new AuxField*[3];
  __int_ke_vec.resize(3);
  for(int ii = 0; ii < 3; ii++) {
    vort[ii]  = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    visc[ii]  = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    rhs[ii]   = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    phi[ii]   = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    phi_c[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    u_o[ii]   = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    __int_ke_vec[ii]   = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    __int_err_prev[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    __visc_rhs_tmp[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
  }
  tmp = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
  __int_ke_tmp   = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
  __int_err_temp = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
  __divg_rhs_tmp = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
  // modal matrix system for the pressure solve;
  // hijack the scalar equation solve for this (with homogenous neumann bcs)
  Femlib::ivalue("N_TIME", 1);
  Femlib::ivalue("N_TIME", ntime);

  alloc_rk = false;
}

void assert_axial_bcs(AuxField** u_arr) {
  int       np     = Geometry::nP();
  int       np2    = np * np;
  real_t*   plane;
  int       plane_j;
  int       node_j;

  AuxField::couple(u_arr[1], u_arr[2], FORWARD);

  for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
    plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;

    // axial velocity
    plane = u_arr[0]->plane(plane_i);
    for(int element_i = 0; element_i < Femlib::ivalue("NELS_X"); element_i++) {
      for(int node_i = 0; node_i < np; node_i++) {
        node_j = element_i * np2 + node_i;

        if(plane_j > 1) plane[node_j] = 0.0;
        if(plane_j ==1) plane[node_j] = 0.0;
      }
    }
    // \tilde{v} velocity
    plane = u_arr[1]->plane(plane_i);
    for(int element_i = 0; element_i < Femlib::ivalue("NELS_X"); element_i++) {
      for(int node_i = 0; node_i < np; node_i++) {
        node_j = element_i * np2 + node_i;

        plane[node_j] = 0.0;
        if(plane_j ==1) plane[node_j] = 0.0;
      }
    }
    // tilde{w} velocity
    plane = u_arr[2]->plane(plane_i);
    for(int element_i = 0; element_i < Femlib::ivalue("NELS_X"); element_i++) {
      for(int node_i = 0; node_i < np; node_i++) {
        node_j = element_i * np2 + node_i;

        if(plane_j > 3) plane[node_j] = 0.0;
        if(plane_j < 2) plane[node_j] = 0.0;
        if(plane_j ==1) plane[node_j] = 0.0;
      }
    }
  }
  // zero out all the nyquist data also
  if(!Geometry::procID()) {
    for(int field_i = 0; field_i < 3; field_i++) {
      plane = u_arr[field_i]->plane(1);
      for(int element_i = 0; element_i < Geometry::nElmt(); element_i++) {
        for(int node_i = 0; node_i < np2; node_i++) {
          plane[element_i*np2+node_i] = 0.0;
        }
      }
    }
  }
  AuxField::couple(u_arr[1], u_arr[2], INVERSE);
for(int field_i = 0; field_i < 3; field_i++) u_arr[field_i]->zeroNyquist();
}

void assert_wall_bcs(AuxField** ui) {
  const int nex  = Femlib::ivalue("NELS_X");
  const int ney  = Femlib::ivalue("NELS_Y");
  const int el_o = nex*(ney - 1);
  const int np   = Geometry::nP();
  const int np2  = np*np;

  for(int field_i = 0; field_i < 3; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      for(int el_i = el_o; el_i < nex*ney; el_i++) {
        for(int pt_i = el_i*np2 + np*(np - 1); pt_i < (el_i+1)*np2; pt_i++) {
          ui[field_i]->plane(plane_i)[pt_i] = 0.0;
        }
      }
    }
  }
}

// compute grad ui, for vector ui
void calc_tensor(AuxField** ui) {
  for(int ii = 0; ii < 3; ii++) {
    for(int jj = 0; jj < 3; jj++) {
      *du[3*ii+jj] = *ui[ii];
      du[3*ii+jj]->gradient(jj);
      if(Geometry::cylindrical() && jj==2) du[3*ii+jj]->divY();
    }
  }
}

// compute the curl of a vector field, ui as wi
void calc_curl(AuxField** ui, AuxField** wi) {
  calc_tensor(ui);

  if(Geometry::cylindrical()) {
    *wi[0]  = *ui[2];
    wi[0]->divY();
  } else {
    *wi[0] = 0.0;
  }
  *wi[0] += *du[2*3+1];
  *wi[0] -= *du[1*3+2];

  *wi[1]  = *du[0*3+2];
  *wi[1] -= *du[2*3+0];

  *wi[2]  = *du[1*3+0];
  *wi[2] -= *du[0*3+1];
}

// radius x divergence (for rhs of pressure poisson equation)
void divg_rhs(Domain* domain, AuxField** ui, AuxField* div, bool mul_rad, bool do_print) {
  double int_div_sq;
  *div = 0.0;

  if(mul_rad) {
    *__divg_rhs_tmp = *ui[0];
//    __divg_rhs_tmp->mulY();
    __divg_rhs_tmp->gradient(0);
    __divg_rhs_tmp->mulY();
    *div += *__divg_rhs_tmp;

    *__divg_rhs_tmp = *ui[1];
    __divg_rhs_tmp->mulY();
    __divg_rhs_tmp->gradient(1);
    *div += *__divg_rhs_tmp;

    *__divg_rhs_tmp = *ui[2];
    __divg_rhs_tmp->gradient(2);
    *div += *__divg_rhs_tmp;
  } else {
    *__divg_rhs_tmp = *ui[0];
    __divg_rhs_tmp->gradient(0);
    *div += *__divg_rhs_tmp;

    *__divg_rhs_tmp = *ui[1];
    __divg_rhs_tmp->mulY();
    __divg_rhs_tmp->gradient(1);
    __divg_rhs_tmp->divY();
    *div += *__divg_rhs_tmp;

    *__divg_rhs_tmp = *ui[2];
    __divg_rhs_tmp->gradient(2);
    __divg_rhs_tmp->divY();
    *div += *__divg_rhs_tmp;
  }

  if(domain) domain->u[0]->smooth(div);

  if(do_print) {
    *__divg_rhs_tmp = *div;
    __divg_rhs_tmp->transform(INVERSE);
    __divg_rhs_tmp->times(*__divg_rhs_tmp, *__divg_rhs_tmp);
    __divg_rhs_tmp->transform(FORWARD);
    int_div_sq = __divg_rhs_tmp->integral();
    if(!Geometry::procID()) cout << "|du.du|^{1/2}: " << sqrt(int_div_sq) << endl;
  }
}

void conv_rhs(Domain* domain, AuxField** ui, AuxField** ni, FieldForce* FF) {
  for(int ii = 0; ii < 3; ii++) ui[ii]->transform(INVERSE);

  for(int ii = 0; ii < 3; ii++) {
    *ni[ii] = 0.0;
    for(int jj = 0; jj < 3; jj++) {
      *__divg_rhs_tmp = *ui[ii];
      if(jj == 2) __divg_rhs_tmp->transform(FORWARD);
      __divg_rhs_tmp->gradient(jj);
      if(jj == 2) __divg_rhs_tmp->divY();
      if(jj == 2) __divg_rhs_tmp->transform(INVERSE);
      if(jj == 2) __divg_rhs_tmp->zeroNyquist();
      *__visc_rhs_tmp[jj] = *ui[jj];
      ni[ii]->timesMinus(*__visc_rhs_tmp[jj], *__divg_rhs_tmp); // -ve of nonlinear term
    }
  }
  *__visc_rhs_tmp[1] = *ui[1];
  *__visc_rhs_tmp[2] = *ui[2];

  __divg_rhs_tmp->times(*__visc_rhs_tmp[2], *__visc_rhs_tmp[2]);
  __divg_rhs_tmp->divY();
  *ni[1] += *__divg_rhs_tmp;                                    // -ve of nonlinear term

  __divg_rhs_tmp->times(*__visc_rhs_tmp[1], *__visc_rhs_tmp[2]);
  __divg_rhs_tmp->divY();
  *ni[2] -= *__divg_rhs_tmp;                                    // -ve of nonlinear term

  for(int ii = 0; ii < 3; ii++) ui[ii]->transform(FORWARD);
  for(int ii = 0; ii < 3; ii++) ui[ii]->zeroNyquist();

  for(int ii = 0; ii < 3; ii++) ni[ii]->transform(FORWARD);
  for(int ii = 0; ii < 3; ii++) ni[ii]->zeroNyquist();

  for(int ii = 0; ii < 3; ii++) *vel0[ii] = *ui[ii];
  *__visc_rhs_tmp[0] = 0.0;
  if(FF) FF->addFourier(__visc_rhs_tmp[0], __divg_rhs_tmp, 0, vel0);
  __visc_rhs_tmp[0]->divY();
  *ni[0] += *__visc_rhs_tmp[0];

  for(int ii = 0; ii < 3; ii++) domain->u[0]->smooth(ni[ii]);
}

void y_conv_rhs(Domain* domain, AuxField** ui, AuxField** ni, FieldForce* FF) {
  for(int ii = 0; ii < 3; ii++) ui[ii]->transform(INVERSE);

  for(int ii = 0; ii < 3; ii++) {
    *ni[ii] = 0.0;
    for(int jj = 0; jj < 3; jj++) {
      *__divg_rhs_tmp = *ui[ii];
      if(jj == 2) __divg_rhs_tmp->transform(FORWARD);
      __divg_rhs_tmp->gradient(jj);
      if(jj == 2) __divg_rhs_tmp->transform(INVERSE);
      if(jj == 2) __divg_rhs_tmp->zeroNyquist();
      if(jj <  2) __divg_rhs_tmp->mulY();
      *__visc_rhs_tmp[jj] = *ui[jj];
      ni[ii]->timesMinus(*__visc_rhs_tmp[jj], *__divg_rhs_tmp); // -ve of nonlinear term
    }
  }
  *__visc_rhs_tmp[1] = *ui[1];
  *__visc_rhs_tmp[2] = *ui[2];

  __divg_rhs_tmp->times(*__visc_rhs_tmp[2], *__visc_rhs_tmp[2]);
  *ni[1] += *__divg_rhs_tmp;                                    // -ve of nonlinear term

  __divg_rhs_tmp->times(*__visc_rhs_tmp[1], *__visc_rhs_tmp[2]);
  *ni[2] -= *__divg_rhs_tmp;                                    // -ve of nonlinear term

  for(int ii = 0; ii < 3; ii++) ui[ii]->transform(FORWARD);
  for(int ii = 0; ii < 3; ii++) ui[ii]->zeroNyquist();

  for(int ii = 0; ii < 3; ii++) ni[ii]->transform(FORWARD);
  for(int ii = 0; ii < 3; ii++) ni[ii]->zeroNyquist();

  for(int ii = 0; ii < 3; ii++) *vel0[ii] = *ui[ii];
  if(FF) FF->addFourier(ni[0], __divg_rhs_tmp, 0, vel0);

  for(int ii = 0; ii < 3; ii++) ni[ii]->divY();

  for(int ii = 0; ii < 3; ii++) domain->u[0]->smooth(ni[ii]);
}

double integrate_ke(AuxField** __u) {
  for(int ii = 0; ii < 3; ii++) {
    *__int_ke_vec[ii] = *__u[ii];
    __int_ke_vec[ii]->transform(INVERSE);
  }
  __int_ke_tmp->innerProduct(__int_ke_vec, __int_ke_vec);
  __int_ke_tmp->transform(FORWARD);
  return __int_ke_tmp->integral();
}

void const_mass_flux_correction(Domain* D, AuxField** ui, real_t c_i) {
  real_t L_x       = Femlib::value("XMAX");
  real_t _refQ     = Femlib::value("Q_BAR");
  real_t getQ, dP;

  if(fabs(_refQ) < 1.0e-6) return;

  if(!Geometry::procID()) {
    getQ  = 2.0 * M_PI * ui[0]->integral(0) / Femlib::ivalue("BETA");
    getQ /= (M_PI * 1.0 * 1.0 * L_x) / Femlib::ivalue("BETA");
    dP    = (_refQ - getQ) / D->Qg;
  }
  MPI_Bcast(&dP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for(int ii = 0; ii < 3; ii++) ui[ii]->axpy(c_i * dP, *D->grn[ii]);
}

void _axpy(AuxField** u1, real_t alpha, AuxField** u2, AuxField** u3) {
  int       np     = Geometry::nP();
  int       np2    = np * np;
  real_t*   plane_1;
  real_t*   plane_2;
  real_t*   plane_3;
  int       plane_j;
  int       node_j;
  int       nex    = Femlib::ivalue("NELS_X");
  int       ney    = Femlib::ivalue("NELS_Y");

  AuxField::couple(u1[1], u1[2], FORWARD);
  AuxField::couple(u2[1], u2[2], FORWARD);
  if(u3 != u1) AuxField::couple(u3[1], u3[2], FORWARD);

  for(int field_i = 0; field_i < 3; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;

      plane_1 = u1[field_i]->plane(plane_i);
      plane_2 = u2[field_i]->plane(plane_i);
      plane_3 = u3[field_i]->plane(plane_i);

      for(int element_i = 0; element_i < Geometry::nElmt(); element_i++) {
        for(int node_i = 0; node_i < np2; node_i++) {
          node_j = element_i * np2 + node_i;

          // nyquist frequency
          if(plane_j == 1) {
            plane_1[node_j] = 0.0;
            continue;
          }
          // wall - homogenous dirichlet bcs
          if(element_i/nex == ney-1 && node_i/np == np-1) {
            plane_1[node_j] = 0.0;
            continue;
          }
          // axis
          if(element_i/nex == 0 && node_i/np == 0) {
            if(plane_j > 3) { // k > 1
              plane_1[node_j] = 0.0;
              continue;
            }
            if(field_i == 0 && plane_j > 1) {
              plane_1[node_j] = 0.0;
              continue;
            }
            if(field_i == 1) {
              plane_1[node_j] = 0.0;
              continue;
            }
            if(field_i == 2 && plane_j < 2) {
              plane_1[node_j] = 0.0;
              continue;
            }
          }
          plane_1[node_j] = alpha * plane_2[node_j] + plane_3[node_j];
        }
      }
    }
  }
  AuxField::couple(u1[1], u1[2], INVERSE);
  AuxField::couple(u2[1], u2[2], INVERSE);
  if(u3 != u1) AuxField::couple(u3[1], u3[2], INVERSE);
}

void zero_mean_radial_pressure_gradient_at_axis(AuxField* dpdy) {
  int       np     = Geometry::nP();
  int       np2    = np * np;
  int       plane_j;
  int       node_j;
  real_t*   plane;
  int       nex    = Femlib::ivalue("NELS_X");
  int       ney    = Femlib::ivalue("NELS_Y");

  for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
    plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
    plane   = dpdy->plane(plane_i);
    // axis
    for(int element_i = 0; element_i < nex; element_i++) {
      for(int node_i = 0; node_i < np; node_i++) {
        node_j = element_i * np2 + node_i;
        if(plane_j < 2) plane[node_j] = 0.0;
      }
    }
    // wall
    for(int element_i = nex*(ney-1); element_i < nex*ney; element_i++) {
      for(int node_i = np*(np-1); node_i < np2; node_i++) {
        node_j = element_i * np2 + node_i;
        plane[node_j] = 0.0;
      }
    }
  }
}

#define ID_DIAGNOSTIC 1

#ifdef ID_DIAGNOSTIC
bool alloc_diagnostics = true;
vector<AuxField*> vcty;
AuxField* enst;
AuxField* pres;

void diagnostics(Domain* domain, bool to_file) {
  double prod, diss, tote, tote_prime, int_dudy, divg;
  Vector dudy;
  ofstream file;

  // allocate if not already done
  if(alloc_diagnostics) {
    alloc_diagnostics = false;
    vel0.resize(3);
    vcty.resize(3);
    __diag_vel_tmp = new AuxField*[3];
    for(int ii = 0; ii < 3; ii++) {
      vel0[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
      *vel0[ii] = *domain->u[ii];
      vel0[ii]->transform(INVERSE);
      vcty[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
      __diag_vel_tmp[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    }
    enst = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    pres = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    return;
  }

  // pressure field was stomped on by the fieldforce, copy back from temporary variable
  pres->transform(INVERSE);
  // transform state into physical space in order to perform pointwise multiplications
  for(int ii = 0; ii < 3; ii++) *__diag_vel_tmp[ii] = *domain->u[ii];
  for(int ii = 0; ii < 3; ii++) __diag_vel_tmp[ii]->transform(INVERSE);

  // energy production, compute as: I = \int_{V} DIV(pu) dV
if(Geometry::cylindrical()) {
  *enst = 0.0;
  for(int ii = 0; ii < 3; ii++) {
    *vcty[ii] = *pres;
    if(ii == 2) vcty[ii]->transform(FORWARD);
    vcty[ii]->gradient(ii);
    if(ii == 2) vcty[ii]->transform(INVERSE);
    if(ii == 2) vcty[ii]->divY();

    if(fabs(Femlib::value("Q_BAR")) > 1.0e-6) {
      // constant mass flux forcing
      if(ii == 0) {
        *enst = *domain->u[PIND]; // use the pressure field to compute the velocity gradient at the wall, then replace with original value
        *domain->u[PIND] = *domain->u[0];
        domain->u[PIND]->transform(FORWARD);
        domain->u[PIND]->gradient(1);
        dudy = Field::normTraction(domain->u[PIND]);
        int_dudy = -2.0 * dudy.y * Femlib::value("KINVIS") / Femlib::value("XMAX");
        MPI_Bcast(&int_dudy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        *vcty[ii] += int_dudy;
        *domain->u[PIND] = *enst;
        *enst = 0.0;
      }
    } else {
      // constant pressure gradient forcing
      if(ii == 0) *vcty[ii] += 4.0*Femlib::value("KINVIS");
    }
    *vcty[ii] *= *__diag_vel_tmp[ii];
    *enst += *vcty[ii];
  }
  enst->transform(FORWARD);
  prod = enst->integral() / Femlib::value("KINVIS");
} else {
  *enst = *domain->u[PIND]; // use the pressure field to compute the velocity gradient at the wall, then replace with original value
  *domain->u[PIND] = *domain->u[0];
  //domain->u[PIND]->transform(FORWARD);
  domain->u[PIND]->gradient(1);
  dudy = Field::normTraction(domain->u[PIND]);
  *domain->u[PIND] = *enst;
  prod = -1.0*dudy.y;// * Femlib::value("XMAX") * Femlib::value("BETA") / Femlib::value("TWOPI");
}

  // energy dissipation, compute as: D = \int_{V} [v^2 + w^2 - 2 w dv/\theta + 2 v dwd/theta]/y^2 + |GRAD u|^2 + |GRAD v|^2 + |GRAD w|^2 dV
  // first do the |GRAD U|^2, U = {u,v,w} terms
  *enst = 0.0;
if(Geometry::cylindrical()) {
  for(int ii = 0; ii < 3; ii++) {
    for(int jj = 0; jj < 3; jj++) {
      if(jj == ii) continue;

      *vcty[0] = *__diag_vel_tmp[ii];
      if(jj == 2) vcty[0]->transform(FORWARD);
      vcty[0]->gradient(jj);
      if(jj == 2) vcty[0]->transform(INVERSE);
      if(jj == 2) vcty[0]->divY();

      *vcty[0] *= *vcty[0];
      *enst += *vcty[0];
    }

    int jj = (ii+1)%3;
    int kk = (ii+2)%3;

    *vcty[0] = *__diag_vel_tmp[jj];
    *vcty[1] = *__diag_vel_tmp[kk];

    if(kk == 2) vcty[0]->transform(FORWARD);
    vcty[0]->gradient(kk);
    if(kk == 2) vcty[0]->transform(INVERSE);
    if(kk == 2) vcty[0]->divY();

    if(jj == 2) vcty[1]->transform(FORWARD);
    vcty[0]->gradient(jj);
    if(jj == 2) vcty[1]->transform(INVERSE);
    if(jj == 2) vcty[1]->divY();

    *vcty[0] *= *vcty[1];
    *vcty[0] *= 2.0;
    *enst -= *vcty[0];
  }

  *vcty[0] = *__diag_vel_tmp[2];
  vcty[0]->gradient(1);
  vcty[0]->divY();
  *vcty[0] *= *__diag_vel_tmp[2];
  *vcty[0] *= 2.0;
  *enst += *vcty[0];

  *vcty[0] = *__diag_vel_tmp[1];
  vcty[0]->transform(FORWARD);
  vcty[0]->gradient(2);
  vcty[0]->transform(INVERSE);
  vcty[0]->divY();
  *vcty[0] *= *__diag_vel_tmp[2];
  vcty[0]->divY();
  *vcty[0] *= 2.0;
  *enst -= *vcty[0];

  // integeate in fourier space
  enst->transform(FORWARD);
  diss = enst->integral();
} else {
  for(int ii = 0; ii < 3; ii++) __diag_vel_tmp[ii]->transform(FORWARD);
  calc_curl(__diag_vel_tmp, vort);
  for(int ii = 0; ii < 3; ii++) __diag_vel_tmp[ii]->transform(INVERSE);
  for(int ii = 0; ii < 3; ii++) vort[ii]->transform(INVERSE);
  for(int ii = 0; ii < 3; ii++) *vcty[ii] = *vort[ii];
  enst->innerProduct(vcty, vcty) *= 0.5;
  enst->transform(FORWARD);
  diss = enst->integral();
}

  // integrate the total energy
  for(int ii = 0; ii < 3; ii++) *vcty[ii] = *__diag_vel_tmp[ii];
  //*vcty[0] += Femlib::value("WAVE_SPEED");
  *enst = 0.0;
  enst->innerProduct(vcty, vcty) *= 0.5;
  enst->transform(FORWARD);
  tote = enst->integral();

  // perturbation kinetic energy
  for(int ii = 0; ii < 3; ii++) {
    *vcty[ii]  = *__diag_vel_tmp[ii];
    *vcty[ii] -= *vel0[ii];
  }
  enst->innerProduct(vcty, vcty) *= 0.5;
  enst->transform(FORWARD);
  tote_prime = enst->integral();

  // divergence
  *enst = 0.0;
  for(int ii = 0; ii < 3; ii++) {
    *vcty[ii] = *__diag_vel_tmp[ii];
    vcty[ii]->transform(FORWARD);
    if(ii == 1) {
      vcty[ii]->divY();
      *enst += *vcty[ii];
      vcty[ii]->mulY();
    }
    vcty[ii]->gradient(ii);
    if(ii == 2) vcty[ii]->divY();
    *enst += *vcty[ii];
  }
  enst->transform(INVERSE);
  *vcty[0] = *enst;
  *enst *= *vcty[0];
  enst->transform(FORWARD);
  divg = enst->integral();
  divg = sqrt(divg);

  // transform state back into fourier space
  pres->transform(FORWARD);
  pres->zeroNyquist();

  if(!Geometry::procID()) {
    if(to_file) {
      file.open("production_dissipation.txt", ios::app);
      file.precision(12);
      file << domain->step << "\t" << prod << "\t" << diss << "\t" << tote << "\t" << tote_prime << "\t" << divg << "\n";
      file.close();
    } else {
      cout << domain->step << "\t" << prod << "\t" << diss << "\t" << tote << "\t" << tote_prime << "\t" << divg << "\n";
    }
  }
}
#endif

double global_min(AuxField* ui, int step) {
  double min = 1.0, g_min;
  const int np   = Geometry::nP();
  const int np2  = np*np;
  real_t* plane;
  double int_u   = ui->integral();
  //double* min_proc = new double[Geometry::nProc()];
  //double min_tmp = 1.0;
  //int min_i;
  int element_i, node_i;

  *__diag_vel_tmp[0] = *ui;
  __diag_vel_tmp[0]->transform(INVERSE);
  for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
    plane = __diag_vel_tmp[0]->plane(plane_i);
    for(int pt_i = 0; pt_i < np2*Geometry::nElmt(); pt_i++) {
      if(plane[pt_i] < min) min = plane[pt_i];
    }
  }
  MPI_Allreduce(&min, &g_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  if(!Geometry::procID()) cout << step << ":\tglobal min: " << g_min << "\tintegral: " << int_u << endl;

/*
  min_proc[Geometry::procID()] = min;
  MPI_Bcast(&min_proc[Geometry::procID()], 1, MPI_DOUBLE, Geometry::procID(), MPI_COMM_WORLD);
  for(int proc_i = 0; proc_i < Geometry::nProc(); proc_i++) {
    if(Geometry::procID()==0) cout << min_proc[proc_i] << endl;
    if(min_proc[proc_i] < min_tmp) {
      min_i   = proc_i;
      min_tmp = min_proc[proc_i];
    }
  }
  if(!Geometry::procID()) cout << step << ":\tglobal min: " << min_tmp << "\tprocessor: " << min_i << endl;
*/

  min = 1.0;
  for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
    plane = __diag_vel_tmp[0]->plane(plane_i);
    for(int pt_i = 0; pt_i < np2*Geometry::nElmt(); pt_i++) {
      if(plane[pt_i] < min) {
        min = plane[pt_i];
        element_i = pt_i/np2;
        node_i = pt_i%np2;
      }
    }
  }
  ROOTONLY cout << Geometry::procID() << "\tmin: " << min 
       << ":\telement x: " << element_i%Femlib::ivalue("NELS_X") 
       << ":\telement y: " << element_i/Femlib::ivalue("NELS_X") 
       << "\tnode x: " << node_i%np 
       << "\tnode y: " << node_i/np 
       << endl;

  //delete[] min_proc;

  return g_min;
}

void zero_nyquist(vector<Field*> ui) {
  const int np   = Geometry::nP();
  const int np2  = np*np;
  real_t*   plane;

  if(!Geometry::procID()) {
    for(int field_i = 0; field_i < 3; field_i++) {
      plane = ui[field_i]->plane(1);
      for(int pt_i = 0; pt_i < np2*Geometry::nElmt(); pt_i++) plane[pt_i] = 0.0;
    }
  }
}

//#define ROTATIONAL_CORRECTION 1

void pc_loop(void (*advection) (Domain*    , 
                                BCmgr*     ,
                                AuxField** , 
                                AuxField** ,
                                FieldForce*),
                Domain*      D ,
                BCmgr*       B ,
                DNSAnalyser* A ,
                FieldForce*  FF, Msys** MMS, AuxField** Us1, AuxField** Uf1, int step)
{
  int           done    = 0;
  const int_t   nmodes  = Geometry::nModeProc();
  const int_t   base    = Geometry::baseMode();
  const real_t  beta    = Femlib::value("BETA");
  const real_t  dt      = Femlib::value("D_T");
  const int_t   ntime   = Femlib::ivalue("N_TIME");
  const real_t  lambda2 = 1.0 / Femlib::value("D_T") / Femlib::value("KINVIS");
  real_t        int_u2, int_du2, int_p, int_dp;
  int           iter    = 0;

  Femlib::ivalue("N_TIME", 1);

  // on first substep put the inital velocity fields in the previous velocity
  // fields storage space for use in the velocity correction scheme once we're done here
  for(int ii = 0; ii < 3; ii++) {
    *Us1[ii] = *D->u[ii];
    if(ii < 2) Us1[ii]->mulY();
  }

  for(int ii = 0; ii < 3; ii++) *u_o[ii] = *D->u[ii];
  for(int ii = 0; ii < 3; ii++) D->u[0]->smooth(u_o[ii]);
  for(int ii = 0; ii < 3; ii++) *u_i[ii] = *D->u[ii];
  for(int ii = 0; ii < 9; ii++) *F_i[ii] = 0.0;
  for(int ii = 0; ii < 9; ii++) F_i[ii]->zeroNyquist();

  // compute the initial ke
  int_u2 = integrate_ke(u_o);
  if(!Geometry::procID()) cout << "initial ke: " << int_u2 << endl;
  int_p = D->u[CIND]->integral();
  if(!Geometry::procID()) cout << "global pressure integral: " << int_p << endl;

  // first step, do a poisson solve for the initial pressure
/*
  if(alloc_pc) {
    divg_rhs(NULL, u_o, rhs[0], true, false); // scale rhs by y
    *rhs[0] /= dt;
    D->u[CIND]->solve(rhs[0], mms);
    divg_rhs(D, u_o, rhs[0], false, false);
    D->u[CIND]->axpy(-1.0*Femlib::value("KINVIS"), *rhs[0]);
    *phi[0] = *D->u[CIND];
  }
*/
  *D->u[CIND] = 0.0;
//*phi[0] = 0.0;
//D->u[CIND]->zeroNyquist();

  while(!done) {
    if(!Geometry::procID()) cout << "iteration: " << iter << endl;
for(int ii = 0; ii < 3; ii++) D->u[0]->smooth(u_o[ii]);
//global_min(D->u[0], 0);
    // 1. viscous step
    //    -- copy the pressure into a temporary buffer as this gets stomped on by the advection routine
    *phi[0] = *D->u[CIND];
    //    -- nonlinear term (previous time level)
//    for(int ii = 0; ii < 3; ii++) *D->u[ii] = *u_o[ii];
//    for(int ii = 0; ii < 3; ii++) *rhs[ii]  = *u_o[ii];
//    advection(D, B, &rhs[0], &F_i[0], FF); // note: axial and radial components are scaled by y here...
    y_conv_rhs(D, u_o, &F_i[0], NULL);
    // on first substep put the inital nonlinear forcing (scaled by y) in the previous nonlinear
    // terms storage space for use in the velocity correction scheme once we're done here
    if(!iter) for(int ii = 0; ii < 3; ii++) *Uf1[ii] = *F_i[ii];
    //    -- nonlinear term (current time level)
//    for(int ii = 0; ii < 3; ii++) *D->u[ii] = *u_i[ii];
//    for(int ii = 0; ii < 3; ii++) *rhs[ii]  = *u_i[ii];
//    advection(D, B, &rhs[0], &F_i[3], FF); // note: axial and radial components are scaled by y here...
    y_conv_rhs(D, u_i, &F_i[3], NULL);
    //    -- copy the pressure back from the buffer
    *D->u[CIND] = *phi[0];
    for(int ii = 0; ii < 3; ii++) {
      //  -- pressure gradient forcing
      *rhs[ii] = *D->u[CIND];
      rhs[ii]->gradient(ii);
      if(ii <  2) rhs[ii]->mulY();
      *rhs[ii] *= (-1.0*dt);
      //  -- previous time step
      *tmp = *u_o[ii];
      tmp->mulY();
      rhs[ii]->axpy(1.0, *tmp);
      //  -- time centered nonlinear term
      if(ii == 2) F_i[0+ii]->mulY();
      if(ii == 2) F_i[3+ii]->mulY();
//      rhs[ii]->axpy(0.5*dt, *F_i[0+ii]);
//      rhs[ii]->axpy(0.5*dt, *F_i[3+ii]);
      rhs[ii]->axpy(dt, *F_i[0+ii]);
      *rhs[ii] *= (-1.0*lambda2);
    }
    // couple the rhs terms so that we are in tilde variables - TODO: viscous term causes blow up...
    for(int ii = 0; ii < 3; ii++) D->u[ii]->evaluateBoundaries(NULL, step, true);
    Field::coupleBCs(D->u[1], D->u[2], FORWARD);
    AuxField::couple(D->u[1], D->u[2], FORWARD);
    AuxField::couple(rhs[1] , rhs[2],  FORWARD);
    for(int ii = 0; ii < 3; ii++) D->u[ii]->solve(rhs[ii], MMS[ii]);
    // uncouple back into hat variables
    AuxField::couple(D->u[1], D->u[2], INVERSE);
//global_min(D->u[0], 1);
//    for(int ii = 0; ii < 3; ii++) D->u[ii]->zeroNyquist();//fails on first time step without this!!
    zero_nyquist(D->u);
//global_min(D->u[0], 2);
/*
for(int ii = 0; ii < 3; ii++) {
*rhs[ii] = *D->u[CIND];
rhs[ii]->gradient(ii);
if(ii == 1) zero_mean_radial_pressure_gradient_at_axis(rhs[ii]);
if(ii == 2) rhs[ii]->divY();
*rhs[ii] *= (-1.0*dt);
rhs[ii]->axpy(1.0, *u_o[ii]);
}

assert_wall_bcs(rhs);
assert_axial_bcs(rhs);
if(!Geometry::procID()) {
for(int el_i = 0; el_i < Femlib::ivalue("NELS_X"); el_i++)
for(int pt_i = el_i*Geometry::nTotElmt(); pt_i < el_i*Geometry::nTotElmt()+Geometry::nP(); pt_i++)
rhs[0]->plane(0)[pt_i] = D->u[0]->plane(0)[pt_i];
}
for(int ii = 0; ii < 3; ii++) rhs[ii]->zeroNyquist();
for(int ii = 0; ii < 3; ii++) *D->u[ii] = *rhs[ii];
*/

    // 2. pressure step
    *phi[0] = *D->u[CIND];
    for(int ii = 0; ii < 3; ii++) *v_i[ii] = *D->u[ii];
    divg_rhs(NULL, v_i, rhs[0], true, false); // scale rhs by y
    *rhs[0] /= dt;
//D->u[0]->smooth(rhs[0]);
    for(int ii = 0; ii < 3; ii++) *v_i[ii] = 0.0;
    B->maintainFourier(1, D->u[PIND], const_cast<const AuxField**>(v_i), const_cast<const AuxField**>(v_i));
    D->u[CIND]->evaluateBoundaries(D->u[CIND], 0);
*rhs[1] = *rhs[0];
    D->u[CIND]->solve(rhs[0], MMS[CIND]);
    D->u[CIND]->zeroNyquist();
    //D->u[0]->smooth(D->u[CIND]);
//*vcty[0] = *D->u[CIND];
//vcty[0]->transform(INVERSE);
//*vcty[1] = 0.0;
//vcty[1]->times(*vcty[0], *vcty[0]);
//vcty[1]->transform(FORWARD);
//int_p = vcty[1]->integral();
//if(!Geometry::procID())cout<<"integral of pressure forcing: "<<int_p<<endl;

    *D->u[CIND] += *phi[0];
    //D->u[0]->smooth(D->u[CIND]);
//global_min(D->u[0], 3);
//*vcty[0] = *D->u[CIND];
//vcty[0]->transform(INVERSE);
//*vcty[1] = 0.0;
//vcty[1]->times(*vcty[0], *vcty[0]);
//vcty[1]->transform(FORWARD);
//int_p = vcty[1]->integral();
//if(!Geometry::procID())cout<<"integral of pressure forcing: "<<int_p<<endl;

    // 3. velocity update
/*
    for(int ii = 0; ii < 3; ii++) {
      *rhs[ii] = *D->u[CIND];
      rhs[ii]->gradient(ii);
      if(ii == 2) rhs[ii]->divY();
      D->u[ii]->axpy(-dt, *rhs[ii]);
      D->u[0]->smooth(D->u[ii]);
    }
*/

    //   -- rotational correction
#ifdef ROTATIONAL_CORRECTION
    for(int ii = 0; ii < 3; ii++) *v_i[ii] = *D->u[ii];
    divg_rhs(D, v_i, rhs[0], false, false);
    D->u[CIND]->axpy(-1.0*Femlib::value("KINVIS"), *rhs[0]);
#endif

    for(int ii = 0; ii < 3; ii++) *v_i[ii] = *D->u[ii];
    const_mass_flux_correction(D, v_i, 1.0);
    for(int ii = 0; ii < 3; ii++) *D->u[ii] = *v_i[ii];

    // 3. convergence test
    int_u2 = integrate_ke(v_i);
    for(int ii = 0; ii < 3; ii++) *v_i[ii] -= *u_i[ii];
    int_du2 = integrate_ke(v_i);
    for(int ii = 0; ii < 3; ii++) *v_i[ii] += *u_i[ii];
    if(!Geometry::procID()) cout << "\t|u|^2: " << int_u2 << "\t|du|^2: " << int_du2 << "\t|du^2/u^2|^{1/2}: " << sqrt(int_du2/int_u2) << "\t";
    int_p = D->u[CIND]->integral();
    if(!Geometry::procID()) cout << "global pressure integral: " << int_p << endl;

    for(int ii = 0; ii < 3; ii++) D->u[ii]->smooth(D->u[ii]);

#ifdef ID_DIAGNOSTIC
    *pres = *D->u[CIND];
    diagnostics(D, true);
#endif
    
    for(int ii = 0; ii < 3; ii++) *u_i[ii] = *D->u[ii];
//global_min(D->u[0], 4);
    //for(int ii = 0; ii < 3; ii++) u_i[ii]->zeroNyquist();

    divg_rhs(D, u_i, rhs[0], false, true); // compute the divergence at the end of the iteration (for testing)

    iter++;
    if(sqrt(int_du2/int_u2) < 1.0e-12) done = 1;
    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

//if(iter==37){
//if(iter==34){
if(iter==400000){
//divg_rhs(D,u_i,rhs[0],false,false);
//*D->u[3]=*rhs[1];
D->dump();abort();
}
  }

/*
  for(int ii = 0; ii < 3; ii++) {
    *rhs[ii] = *D->u[CIND];
    rhs[ii]->gradient(ii);
    if(ii == 2) rhs[ii]->divY();
  }
  _axpy(u_i, -dt, rhs, u_i);
  for(int ii = 0; ii < 3; ii++) D->u[0]->smooth(u_i[ii]);
  assert_axial_bcs(u_i);
  assert_wall_bcs(u_i);
  for(int ii = 0; ii < 3; ii++) *D->u[ii] = *u_i[ii];
*/
  for(int ii = 0; ii < 3; ii++) {
    *rhs[ii] = *D->u[CIND];
    rhs[ii]->gradient(ii);
    if(ii == 2) rhs[ii]->divY();
    D->u[ii]->axpy(-dt, *rhs[ii]);
    D->u[ii]->smooth(D->u[ii]);
  }
{
int_u2 = integrate_ke(v_i);
if(!Geometry::procID()) cout << "end of step, final ke: " << int_u2 << endl;
int_p = D->u[CIND]->integral();
if(!Geometry::procID()) cout << "global pressure integral: " << int_p << endl;
}

  Femlib::ivalue("N_TIME", ntime);

  alloc_pc = false;
}

void pc_loop_no_forcing(
                Domain*      D ,
                BCmgr*       B ,
                DNSAnalyser* A ,
                Msys** MMS, AuxField** Us1, AuxField** Uf1, int step)
{
  int           done    = 0;
  const real_t  dt      = Femlib::value("D_T");
  real_t        int_u2, int_du2, int_p, int_dp;
  int           iter    = 0;

  // on first substep put the inital velocity fields in the previous velocity
  // fields storage space for use in the velocity correction scheme once we're done here
  for(int ii = 0; ii < 3; ii++) {
    *Us1[ii] = *D->u[ii];
    if(ii < 2) Us1[ii]->mulY();
  }

  for(int ii = 0; ii < 3; ii++) *u_o[ii] = *D->u[ii];
  for(int ii = 0; ii < 3; ii++) *u_i[ii] = *D->u[ii];

  if(!Geometry::procID()) cout << "\ndoing pressure correction iterative solve, tracer index: " << CIND << ", pressure index: " << PIND << "\n" << endl;

  // compute the initial ke
  int_u2 = integrate_ke(u_o);
  if(!Geometry::procID()) cout << "initial ke: " << int_u2 << endl;
  int_p = D->u[CIND]->integral();
  if(!Geometry::procID()) cout << "global pressure integral: " << int_p << endl;

  // first step, do a poisson solve for the initial pressure
/*
  if(alloc_pc) {
    divg_rhs(NULL, u_o, rhs[0], true, false); // scale rhs by y
    *rhs[0] /= dt;
    D->u[CIND]->solve(rhs[0], mms);
    divg_rhs(D, u_o, rhs[0], false, false);
    D->u[CIND]->axpy(-1.0*Femlib::value("KINVIS"), *rhs[0]);
    *phi[0] = *D->u[CIND];
  }
*/
//*D->u[CIND] = 0.0;
//*phi[0] = 0.0;

  while(!done) {
    if(!Geometry::procID()) cout << "iteration: " << iter << "\t";
    // 1. viscous step
    // on first substep put the inital nonlinear forcing (scaled by y) in the previous nonlinear
    // terms storage space for use in the velocity correction scheme once we're done here
    for(int ii = 0; ii < 3; ii++) {
      //  -- pressure gradient forcing
      *rhs[ii] = *D->u[CIND];
      rhs[ii]->gradient(ii);
      if(ii == 1) zero_mean_radial_pressure_gradient_at_axis(rhs[ii]);
      if(ii == 2) rhs[ii]->divY();
      *v_i[ii] = *u_o[ii];
      v_i[ii]->axpy(-1.0*dt, *rhs[ii]);
    }
/*
for(int ii = 0; ii < 3; ii++) {
*rhs[ii] = *D->u[CIND];
rhs[ii]->gradient(ii);
if(ii == 1) zero_mean_radial_pressure_gradient_at_axis(rhs[ii]);
if(ii == 2) rhs[ii]->divY();
*rhs[ii] *= (-1.0*dt);
rhs[ii]->axpy(1.0, *u_o[ii]);
}

assert_wall_bcs(rhs);
assert_axial_bcs(rhs);
if(!Geometry::procID()) {
for(int el_i = 0; el_i < Femlib::ivalue("NELS_X"); el_i++)
for(int pt_i = el_i*Geometry::nTotElmt(); pt_i < el_i*Geometry::nTotElmt()+Geometry::nP(); pt_i++)
rhs[0]->plane(0)[pt_i] = D->u[0]->plane(0)[pt_i];
}
for(int ii = 0; ii < 3; ii++) rhs[ii]->zeroNyquist();
for(int ii = 0; ii < 3; ii++) *D->u[ii] = *rhs[ii];
*/

    // 2. pressure step
    *phi[0] = *D->u[CIND];
    divg_rhs(NULL, v_i, rhs[0], true, false); // scale rhs by y
    *rhs[0] /= dt;
    for(int ii = 0; ii < 3; ii++) *v_i[ii] = 0.0;
    B->maintainFourier(1, D->u[PIND], const_cast<const AuxField**>(v_i), const_cast<const AuxField**>(v_i));
    D->u[CIND]->evaluateBoundaries(D->u[CIND], 0);
    //D->u[CIND]->solve(rhs[0], mms);
    D->u[CIND]->solve(rhs[0], MMS[CIND]);
    D->u[CIND]->zeroNyquist();
    *D->u[CIND] += *phi[0];

    // 3. convergence test
    for(int ii = 0; ii < 3; ii++) *v_i[ii] = *D->u[ii];
    int_u2 = integrate_ke(v_i);
    for(int ii = 0; ii < 3; ii++) *v_i[ii] -= *u_i[ii];
    int_du2 = integrate_ke(v_i);
    for(int ii = 0; ii < 3; ii++) *v_i[ii] += *u_i[ii];
    if(!Geometry::procID()) cout << "\t|u|^2: " << int_u2 << "\t|du|^2: " << int_du2 << "\t|du^2/u^2|^{1/2}: " << sqrt(int_du2/int_u2) << "\t";
    int_p = D->u[CIND]->integral();
    if(!Geometry::procID()) cout << "global pressure integral: " << int_p << "\t";

    for(int ii = 0; ii < 3; ii++) D->u[ii]->smooth(D->u[ii]);

#ifdef ID_DIAGNOSTIC
    *pres = *D->u[CIND];
    diagnostics(D, true);
#endif
    
    for(int ii = 0; ii < 3; ii++) *u_i[ii] = *D->u[ii];

    divg_rhs(D, u_i, rhs[0], false, true); // compute the divergence at the end of the iteration (for testing)

    iter++;
    if(sqrt(int_du2/int_u2) < 1.0e-14) done = 1;
    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for(int ii = 0; ii < 3; ii++) {
    *rhs[ii] = *phi_c[0];
    rhs[ii]->gradient(ii);
    if(ii == 2) rhs[ii]->divY();
  }
  _axpy(u_i, -dt, rhs, u_i);
  for(int ii = 0; ii < 3; ii++) D->u[0]->smooth(u_i[ii]);
  for(int ii = 0; ii < 3; ii++) *D->u[ii] = *u_i[ii];

{
int_u2 = integrate_ke(v_i);
if(!Geometry::procID()) cout << "end of step, final ke: " << int_u2 << endl;
int_p = D->u[CIND]->integral();
if(!Geometry::procID()) cout << "global pressure integral: " << int_p << endl;
}

  *D->u[PIND] = *D->u[CIND];

  alloc_pc = false;
}

void integrate (void (*advection) (Domain*    , 
                                   BCmgr*     ,
                                   AuxField** , 
                                   AuxField** ,
                                   FieldForce*),
                Domain*      D ,
                BCmgr*       B ,
                DNSAnalyser* A ,
                FieldForce*  FF)
// ---------------------------------------------------------------------------
// On entry, D contains storage (in the following order!) for:
// -- velocity Fields 'u', 'v' (and 'w' if 2D3C or 3D),
// -- optional scalar Field 'c',
// -- constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NCOM = D -> nVelCmpt();              // -- Number of velocity components.
  //NADV = D -> nAdvect();               // -- Number of advected fields.
  NADV = NCOM;
  NDIM = Geometry::nDim();	       // -- Number of space dimensions.
  NORD = Femlib::ivalue ("N_TIME");    // -- Time integration order.
  C3D  = Geometry::cylindrical() && NDIM == 3;
  PIND = (D->hasScalar()) ? NCOM+1 : NCOM;

  //else if(Femlib::ivalue("RK_SUBSTEP")) { CIND = D->nVelCmpt(); PIND = CIND; }
CIND = D->nVelCmpt();
  if(!Geometry::procID()) cout << "pressure index: " << PIND << endl;
  if(!Geometry::procID()) cout << "  scalar index: " << CIND << endl;

  int_t              i, j, k;
  const real_t       dt    = Femlib:: value ("D_T");
  const int_t        nStep = Femlib::ivalue ("N_STEP");
  const int_t        nZ    = Geometry::nZProc();
  static Msys**      MMS;
  static AuxField*** Us;
  static AuxField*** Uf;
  Field*             Pressure = D -> u[PIND];
  bool substep = true;

  if (!MMS) {			// -- Initialise static storage.

    // -- Create multi-level storage for velocities and forcing.

    const int_t ntot  = Geometry::nTotProc();
    real_t*     alloc = new real_t [static_cast<size_t>(2 * NADV*NORD * ntot)];
    Us                = new AuxField** [static_cast<size_t>(2 * NORD)];
    Uf                = Us + NORD;

    for (k = 0, i = 0; i < NORD; i++) {
      Us[i] = new AuxField* [static_cast<size_t>(2 * NADV)];
      Uf[i] = Us[i] + NADV;
      for (j = 0; j < NADV; j++) {
        Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt);
        Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt);
      }
    }

    // -- Create global matrix systems.

    MMS = preSolve (D);

    // -- Create multi-level storage for pressure BCS.

    B -> buildComputedBCs (Pressure);

    // -- Apply coupling to radial & azimuthal velocity BCs.

    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);
  }

  // -- Because we may restart from scratch on each call, zero these:

  *Pressure = 0.0;

  for (i = 0; i < NORD; i++)
    for (j = 0; j < NADV; j++) {
      *Us[i][j] = 0.0;
      *Uf[i][j] = 0.0;
    }

  // -- Solve the Stokes flow problem with unit forcing
  if(!D->grn[0] && fabs(Femlib::value("Q_BAR")) > 1.0e-6) {
    real_t            L_x   = Femlib::value("XMAX");
    vector<AuxField*> mff_tmp;
    mff_tmp.resize(3);

    for (i = 0; i < 3; i++) {
      D->grn[i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), D->elmt, 'g'+i);
      mff_tmp[i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), D->elmt, 'k'+i);
      *mff_tmp[i] = *D->u[i];
      *D->u[i] = 0.0;
    }

    // -- Set the constant forcing
    for (i = 0; i < NCOM; i++) *Uf[0][i] = 0.0;
    ROOTONLY {
      Uf[0][0] -> addToPlane (0, 1.0);
      if (Geometry::cylindrical()) Uf[0][0] -> mulY ();
    }

    D->step += 1;

    // -- Update high-order pressure BC storage.
    B -> maintainFourier (D -> step, Pressure, const_cast<const AuxField**>(Us[0]), const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (Pressure, D -> step);

    // -- Complete unconstrained advective substep and compute pressure.
    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }

    waveProp (D, const_cast<const AuxField***>(Us), const_cast<const AuxField***>(Uf));
    for (i = 0; i < NADV; i++) AuxField::swapData (D -> u[i], Us[0][i]);

    rollm     (Uf, NORD, NADV);
    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, PIND,  Uf[0][0], MMS[PIND]);

    // -- Correct velocities for pressure.
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.
    for (i = 0; i < NADV; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NADV);

    // -- Re-evaluate velocity (possibly time-dependent) BCs.
    for (i = 0; i < NADV; i++)  {
      D -> u[i] -> evaluateBoundaries (NULL,     D -> step, false);
      D -> u[i] -> bTransform         (FORWARD);
      D -> u[i] -> evaluateBoundaries (Pressure, D -> step, true);
    }
    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

    // -- Viscous correction substep.
    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i < NADV; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if (C3D) AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    D->step -= 1;

    // this will be broadcast to the other procs when applied
    D->Qg  = 2.0 * M_PI * D->u[0]->integral(0) / Femlib::ivalue("BETA");
    D->Qg /= (M_PI * 1.0 * 1.0 * L_x / Femlib::ivalue("BETA"));
    if(!Geometry::procID()) cout << "Stokes + unit forcing volumetric flux: " << D->Qg << endl;
    if(!Geometry::procID()) cout << "                          ux integral: " << 2.0 * M_PI * D->u[0]->integral(0) / Femlib::ivalue("BETA") << endl;
    if(!Geometry::procID()) cout << "                          pipe length: " << L_x << endl;

    // -- Resetting fields
    for (i = 0; i < 3; i++) {
      *D->grn[i] = *D->u[i];
      *D->u[i]   = *mff_tmp[i];
    }

    *Pressure = 0.0;
    for (i = 0; i < NORD; i++)
      for (j = 0; j < NADV; j++) {
        *Us[i][j] = 0.0;
        *Uf[i][j] = 0.0;
      }
    // do we really need to do this again??
    //B -> buildComputedBCs (Pressure);
  }

  // -- The following timestepping loop implements equations (15--18) in [5].

#ifdef ID_DIAGNOSTIC
  // setup only
  diagnostics(D, true);
#endif
  init_fields(D);

  while (D -> step < nStep) {

#ifdef ID_DIAGNOSTIC
    diagnostics(D, true);
#endif

    // iterate on the pressure correction scheme for the first step
    if(substep && Femlib::ivalue("PC_STEP")) {
      for(int ii = 0; ii < Femlib::ivalue("PC_STEP"); ii++) {
        if(!Geometry::procID()) cout << "doing pressure correction, step: " << ii << endl;
        pc_loop(advection, D, B, A, FF, MMS, Us[1], Uf[1], ii);
#ifdef ID_DIAGNOSTIC
        *pres = *D->u[CIND];
        diagnostics(D, true);
#endif
        D->time += dt;
      }
      D->step = Femlib::ivalue("PC_STEP");

      substep = false;

      continue;
    }

    // -- Compute nonlinear terms from previous velocity field.
    //    Add physical space forcing, again at old time level.
    //
    advection (D, B, Us[0], Uf[0], FF);
    
    // -- Now update the time (remainder including BCs at new time level).

    D -> step += 1;
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Update high-order pressure BC storage.

    B -> maintainFourier (D -> step, Pressure,
			  const_cast<const AuxField**>(Us[0]),
			  const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (Pressure, D -> step);

    // -- Complete unconstrained advective substep and compute pressure.

    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }

    waveProp (D, const_cast<const AuxField***>(Us),
	         const_cast<const AuxField***>(Uf));
    for (i = 0; i < NADV; i++) AuxField::swapData (D -> u[i], Us[0][i]);

    rollm     (Uf, NORD, NADV);
    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, PIND,  Uf[0][0], MMS[PIND]);

#ifdef ID_DIAGNOSTIC
    // copy over before this gets stomped on by the fieldforce
    *pres = *D->u[PIND];
#endif

    // -- Correct velocities for pressure.

    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NADV; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NADV);

    // -- Re-evaluate velocity (possibly time-dependent) BCs.

    for (i = 0; i < NADV; i++)  {
      D -> u[i] -> evaluateBoundaries (NULL,     D -> step, false);
      D -> u[i] -> bTransform         (FORWARD);
      D -> u[i] -> evaluateBoundaries (Pressure, D -> step, true);
    }
    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }

    if(D->step == 1) Femlib::ivalue("N_TIME", 1);
    if(D->step == 2) Femlib::ivalue("N_TIME", 2);
    for (i = 0; i < NADV; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if(D->step == 1) Femlib::ivalue("N_TIME", NORD);
    if(D->step == 2) Femlib::ivalue("N_TIME", NORD);
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // constant flow rate
    if(fabs(Femlib::value("Q_BAR")) > 1.0e-6) {
      real_t L_x       = Femlib::value("XMAX");
      real_t _refQ     = Femlib::value("Q_BAR");
      real_t getQ, dP;
      if(!Geometry::procID()) {
        getQ  = 2.0 * M_PI * D->u[0]->integral(0) / Femlib::ivalue("BETA");
        getQ /= (M_PI * 1.0 * 1.0 * L_x) / Femlib::ivalue("BETA");
        dP    = (_refQ - getQ) / D->Qg;
      }
      MPI_Bcast(&dP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      for(i = 0; i < NADV; i++) D->u[i]->axpy(dP, *D->grn[i]);
    }

    // -- Process results of this step.

    A -> analyse (Us[0], Uf[0]);
  }
}

static void waveProp (Domain*           D ,
		      const AuxField*** Us,
		      const AuxField*** Uf)
// ---------------------------------------------------------------------------
// Compute the first substep of stiffly-stable timestepping scheme.
//
// On entry, the most recent velocity fields are in Us, and the most
// recent nonlinear terms in Uf.  The intermediate velocity field u^ is
// computed and left in D's velocity areas.
//
// This is the only routine that makes explicit use of the multi time
// level structure of Us & Uf.
// ---------------------------------------------------------------------------
{
  int_t             i, q;
  vector<AuxField*> H (NADV);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NADV; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const int_t    Je = min (D -> step, NORD);
  vector<real_t> alpha (Integration::OrderMax + 1);
  vector<real_t> beta  (Integration::OrderMax);

  Integration::StifflyStable (Je, &alpha[0]);
  Integration::Extrapolation (Je, &beta [0]);
  Blas::scal (Je, Femlib::value ("D_T"), &beta[0],  1);

  for (i = 0; i < NADV; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[q][i]);
      H[i] -> axpy ( beta [q]    , *Uf[q][i]);
    }
}


static void setPForce (const AuxField** Us,
		       AuxField**       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in Us.  Create div u^ / D_T
// in the first dimension of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t dt = Femlib::value ("D_T");

  for (i = 0; i < NDIM; i++) (*Uf[i] = *Us[i]) . gradient (i);

  for (i = 1; i < NDIM; i++) *Uf[0] += *Uf[i];

  *Uf[0] /= dt;
}


static void project (const Domain* D ,
		     AuxField**    Us,
		     AuxField**    Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in Us.  Constrain velocity field:
//
//                    u^^ = u^ - D_T * grad P,
//
// then scale by -1.0 / (D_T * KINVIS) to create forcing for viscous step
// (this is -1.0 / (D_T  * diffusivity) in the case of a scalar field).
//
// u^^ is left in Uf.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  const real_t beta  =  1.0 / Femlib::value ("KINVIS");
  const real_t Pr    =        Femlib::value ("PRANDTL");
  real_t mass_flux;

  for (i = 0; i < NADV; i++) {
    Field::swapData (Us[i], Uf[i]);
    if (Geometry::cylindrical() && i >= 2) Uf[i] -> mulY();
    *Uf[i] *= alpha;
  }

  // -- For scalar, use diffusivity instead of viscosity.
  if (NADV > NCOM) *Uf[NCOM] *= Pr;

  for (i = 0; i < NDIM; i++) {
    (*Us[0] = *D -> u[PIND]) . gradient (i);

    if (Geometry::cylindrical() && i <  2) Us[0] -> mulY();
    Uf[i] -> axpy (beta, *Us[0]);
  }
}


static Msys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for each Field of D.  If iterative solution
// is selected for any Field, the corresponding ModalMatrixSystem pointer
// is set to zero.
//
// ITERATIVE >= 1 selects iterative solver for velocity components,
// ITERATIVE >= 2 selects iterative solver for non-zero pressure Fourier modes.
// ---------------------------------------------------------------------------
{
  const int_t             nmodes = Geometry::nModeProc();
  const int_t             base   = Geometry::baseMode(); 
  const int_t             itLev  = Femlib::ivalue ("ITERATIVE");
  const real_t            beta   = Femlib:: value ("BETA");
  const vector<Element*>& E = D -> elmt;
  Msys**                  M = new Msys* [static_cast<size_t>(NADV + 1)];
  int_t                   i;

  vector<real_t> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  //real_t   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
real_t   lambda2 = 1.0 / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.
  for (i = 0; i < NCOM; i++) {
    M[i] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[i], (itLev) ? JACPCG : DIRECT);
  }

  // -- Scalar system.
  if (NADV != NCOM) {
    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL");
    M[NCOM] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[NCOM],(itLev < 1)?DIRECT:JACPCG);
  }

  // -- Pressure system.
  if (itLev > 1)
    M[PIND] = new Msys
      (0.0, beta, base, nmodes, E, D -> b[PIND], MIXED);
  else
    M[PIND] = new Msys
      (0.0, beta, base, nmodes, E, D -> b[PIND], DIRECT);

  return M;
}


static void Solve (Domain*     D,
		   const int_t i,
		   AuxField*   F,
		   Msys*       M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for D->u[i], using F as a forcing Field.
// Iterative or direct solver selected on basis of field type, step,
// time order and command-line arguments.
// ---------------------------------------------------------------------------
{
  const int_t step = D -> step;

  if (i < NADV && step < NORD) { // -- We need a temporary matrix system.
    const int_t Je     = min (step, NORD);
    const int_t base   = Geometry::baseMode();
    const int_t nmodes = Geometry::nModeProc();

    vector<real_t> alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t   lambda2 = (i == NCOM) ? // -- True for scalar diffusion.
      alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL") :
      alpha[0] / Femlib::value ("D_T * KINVIS");
    const real_t   beta    = Femlib::value ("BETA");

    Msys* sys_tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, sys_tmp);
    delete sys_tmp;

  } else D -> u[i] -> solve (F, M);
}
