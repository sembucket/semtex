//////////////////////////////////////////////////////////////////////////////
// helmholtz.C:  routines to solve elliptic problems in one variable.
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>


void Helmholtz (Domain*     D      ,
		const char* forcing)
// ---------------------------------------------------------------------------
// Solve Helmholtz's equation
//                                  2
//               div grad u - LAMBDA  u = f(x,y)
//
// subject to BCs.
// ---------------------------------------------------------------------------
{
  const real    lambda2 = Femlib::value ("LAMBDA2");
  const real    beta    = Femlib::value ("BETA");
  const integer nmodes  = Geometry::nModeProc();
  const integer base    = Geometry::baseMode();
  const integer nz      = Geometry::nZProc();
  real*         alloc   = new real [(size_t) Geometry::nTotProc() * nz];
  AuxField*     Force   = new AuxField (alloc, nz, D -> elmt);
  SolverKind    method  = ((int) Femlib::value("ITERATIVE")) ? JACPCG : DIRECT;

  if   (forcing) (*Force = forcing) . transform (+1);
  else            *Force = 0.0;
  
  ModalMatrixSys* M = new ModalMatrixSys
    (lambda2, beta, base, nmodes, D -> elmt, D -> b[0], method);

  D -> u[0] -> solve (Force, M);

  D -> step++;
}
