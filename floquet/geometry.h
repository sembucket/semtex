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
// Geometry has been extended to include the processor ID, number of
// processors, number of data planes per processor, etc.
//
// This version adds a number of functions that are used in
// eigensystem analysis.
//
// $Id$
// ===========================================================================
{
public:
  enum CoordSys { Cartesian, Cylindrical };
  enum Category { O2_2D, O2_3D, O2_3D_SYMM, SO2_3D }; // Refer to README file.

  static void set (const int_t, const int_t);

  static CoordSys    system()       { return _csys;               }
  static Category    problem()      { return _cat;                }
  static const char* symmetry();    // -- Return string corresponding to _cat.
  static bool        cylindrical () { return _csys == Geometry::Cylindrical; }

  static int_t nP        () { return _np;                   }
  static int_t nZ        () { return _nz;                   }
  static int_t nElmt     () { return _nel;                  }
  static int_t nTotElmt  () { return _np * _np;             }
  static int_t nExtElmt  () { return 4 * (_np - 1);         }
  static int_t nIntElmt  () { return (_np - 2) * (_np - 2); }
  static int_t nMode     () { return (_nz + 1) >> 1;        }
  static int_t nDim      () { return _npert;                }
  static int_t kFund     () { return _kfund;                }
  static int_t nPlane    () { return _nel * nTotElmt();     }
  static int_t nBnode    () { return _nel * nExtElmt();     }
  static int_t nInode    () { return _nel * nIntElmt();     }
  static int_t nTot      () { return _nz  * nPlane();       }
  static int_t planeSize () { return _psize;                }
  static int_t nTotal    () { return _nz * _psize;          }

  static int_t nProc     () { return _nproc;                }
  static int_t procID    () { return _pid;                  }
  static int_t nZProc    () { return _nzp;                  }
  static int_t nZ32      () { return (_nproc > 1) ? _nzp : (3 * _nz) >> 1; }
  static int_t nTotProc  () { return _nzp * _psize;         }
  static int_t nModeProc () { return nMode() / _nproc;      }
  static int_t baseMode  () { return _pid * nModeProc();    }
  static int_t basePlane () { return _pid * _nzp;           }
  static int_t nBlock    () { return _psize / _nproc;       }
  
  // -- These are specific to eigensystem analysis:

  static int_t nBase     () { return _nbase;                }
  static int_t nPert     () { return _npert;                }
  static int_t nSlice    () { return _nslice;               }

private:
  static int_t    _nproc ;     // Number of processors.
  static int_t    _pid   ;     // ID for this processor, starting at 0.
  static int_t    _ndim  ;     // Number of space dimensions
  static int_t    _np    ;     // Number of points along element edge.
  static int_t    _nz    ;     // Number of planes (total).
  static int_t    _nzp   ;     // Number of planes per processor.
  static int_t    _nel   ;     // Number of elements.
  static int_t    _psize ;     // nPlane rounded up to suit restrictions.
  static int_t    _kfund ;     // Wavenumber of first non-zero Fourier mode.
  static CoordSys _csys  ;     // Coordinate system (Cartesian/cylindrical).
  static Category _cat   ;     // Problem category.

  // -- These are specific to eigensystem analysis:

  static int_t    _npert ;     // Number of perturbation velocity components.
  static int_t    _nbase ;     // Number of base velocity field components.
  static int_t    _nslice;     // Number of base velocity fields.

};
#endif
