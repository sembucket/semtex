#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cfemdef.h>

class Geometry
// ===========================================================================
// Details of geometric representation used for scalar fields.  Static
// functions make information globally accessible.
//
// NB: this version modified to suit rertictions of dual: 1 process, 3
// planes, 2 modes.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// In all cases, 2D quad elements are employed, with a possible
// extension by Fourier expansions in the third dimension.  While the
// representation implied by this class is not necessarily conforming,
// the same order of interpolation is used in each element.
// Equal-order interpolation is used in each direction on faces of
// quads.
//
// $Id$
// ===========================================================================
{
public:
  enum CoordSys { Cartesian, Cylindrical };

  static void set (const int_t, const int_t, const int_t, const CoordSys);

  static CoordSys system      () { return csys; }  
  static bool     cylindrical () { return csys == Geometry::Cylindrical; }

  static int_t  nP        () { return np;                  }
  static int_t  nZ        () { return nz;                  }
  static int_t  nElmt     () { return nel;                 }
  static int_t  nTotElmt  () { return np * np;             }
  static int_t  nExtElmt  () { return 4 * (np - 1);        }
  static int_t  nIntElmt  () { return (np - 2) * (np - 2); }
  static int_t  nMode     () { return 2;                   }
  static int_t  nDim      () { return 3;                   }
  static int_t  nPlane    () { return nel * nTotElmt();    }
  static int_t  nBnode    () { return nel * nExtElmt();    }
  static int_t  nInode    () { return nel * nIntElmt();    }
  static int_t  nTot      () { return nz  * nPlane();      }
  static int_t  planeSize () { return psize;               }
  static int_t  nTotal    () { return nz * psize;          }

  static int_t  nProc     () { return nproc;               }
  static int_t  procID    () { return pid;                 }
  static int_t  nZProc    () { return nz;                  }
  static int_t  nZ32      () { return (nproc > 1) ? nzp : (3 * nz) >> 1; }
  static int_t  nTotProc  () { return nzp * psize;         }
  static int_t  nModeProc () { return 2;                   }
  static int_t  baseMode  () { return pid * nModeProc();   }
  static int_t  basePlane () { return pid * nzp;           }
  static int_t  nBlock    () { return psize / nproc;       }

private:
  static int_t    nproc ;	// Number of processors.
  static int_t    pid   ;	// ID for this processor, starting at 0.
  static int_t    ndim  ;       // Number of space dimensions.
  static int_t    np    ;	// Number of points along element edge.
  static int_t    nz    ;	// Number of planes (total).
  static int_t    nzp   ;	// Number of planes per processor.
  static int_t    nel   ;	// Number of elements.
  static int_t    psize ;	// nPlane rounded up to suit restrictions.
  static CoordSys csys  ;	// Coordinate system (Cartesian/cylindrical).

};
#endif
