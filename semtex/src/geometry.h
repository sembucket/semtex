#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cfemdef.h>

class Geometry
// ===========================================================================
// Details of geometric representation used for scalar fields.  Static
// functions make information globally accessible.
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
// With the introduction of concurrent execution, the concept of
// geometry has been extended to include the processor ID, number of
// processors, number of data planes per processor, etc.
//
// $Id$
// ===========================================================================
{
public:
  enum CoordSys { Cartesian, Cylindrical };

  static void set (const integer,const integer,const integer,const CoordSys);

  static CoordSys system    () { return _csys;                 }  
  static bool     cylindrical () { return _csys == Geometry::Cylindrical; }

  static integer  nP        () { return _np;                   }
  static integer  nZ        () { return _nz;                   }
  static integer  nElmt     () { return _nel;                  }
  static integer  nTotElmt  () { return _np * _np;             }
  static integer  nExtElmt  () { return 4 * (_np - 1);         }
  static integer  nIntElmt  () { return (_np - 2) * (_np - 2); }
  static integer  nMode     () { return (_nz + 1) >> 1;        }
  static integer  kFund     () { return _kfund;                }
  static integer  nDim      () { return _ndim;                 }
  static integer  nPlane    () { return _nel * nTotElmt();     }
  static integer  nBnode    () { return _nel * nExtElmt();     }
  static integer  nInode    () { return _nel * nIntElmt();     }
  static integer  nTot      () { return _nz  * nPlane();       }
  static integer  planeSize () { return _psize;                }
  static integer  nTotal    () { return _nz * _psize;          }

  static integer  nProc     () { return _nproc;                }
  static integer  procID    () { return _pid;                  }
  static integer  nZProc    () { return _nzp;                  }
  static integer  nZ32      () { return (_nproc > 1) ? _nzp : (3 * _nz) >> 1; }
  static integer  nTotProc  () { return _nzp * _psize;         }
  static integer  nModeProc () { return nMode() / _nproc;      }
  static integer  baseMode  () { return _pid * nModeProc();    }
  static integer  basePlane () { return _pid * _nzp;           }
  static integer  nBlock    () { return _psize / _nproc;       }

private:
  static integer  _nproc ;	// Number of processors.
  static integer  _pid   ;	// ID for this processor, starting at 0.
  static integer  _ndim  ;	// Number of space dimensions.
  static integer  _np    ;	// Number of points along element edge.
  static integer  _nz    ;	// Number of planes (total).
  static integer  _nzp   ;	// Number of planes per processor.
  static integer  _nel   ;	// Number of elements.
  static integer  _psize ;	// nPlane rounded up to suit restrictions.
  static integer  _kfund ;	// Wavenumber of first non-zero Fourier mode.
  static CoordSys _csys  ;	// Coordinate system (Cartesian/cylindrical).

};
#endif
