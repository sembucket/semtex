//////////////////////////////////////////////////////////////////////////////
// helmholtz.C:  routines to solve elliptic problems in one variable.
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <Sem.h>


void Helmholtz (Domain*   D,
		AuxField* F)
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
  SolverKind    method  = ((int) Femlib::value("ITERATIVE")) ? JACPCG : DIRECT;

  ModalMatrixSys* M = new ModalMatrixSys
    (lambda2, beta, base, nmodes, D -> elmt, D -> b[0], method);

  D -> u[0] -> solve (F, M);

  D -> step++;
}
