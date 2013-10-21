//////////////////////////////////////////////////////////////////////////////
// helmholtz.C:  routines to solve elliptic problems in one variable.
//
// Copyright (c) 1994<-->$Date$, Hugh Blackburn
//
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
// 02110-1301 USA
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>


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
  SolverKind    method  = (Femlib::ivalue("ITERATIVE")) ? JACPCG : DIRECT;

  ModalMatrixSys* M = new ModalMatrixSys
    (lambda2, beta, base, nmodes, D -> elmt, D -> b[0], method);

  D -> u[0] -> solve (F, M);

  D -> step++;
}
