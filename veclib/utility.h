#ifndef UTILITY_H
#define UTILITY_H
///////////////////////////////////////////////////////////////////////////////
// Utility.h:  C++ wrapper for veclib utility and memory routines.
//
// Copyright (c) 1994,2003 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <cfemdef>

template<class T> inline T sqr(T x)             { return x * x;            }
template<class T> inline T sgn(T x)             { return (x < 0) ? -1 : 1; }
template<class T> inline T clamp(T t, T a, T b) { return max(min(t,b),a);  }

// Max & Min are part of STL now, so have been removed from the above.

#ifndef M_PI
const double M_PI   = 3.14159265358979323846;
#endif

#ifndef TWOPI
const double TWOPI  = 6.28318530717958647692;
#endif

#ifndef PI_180
const double PI_180 = 0.01745329251994329576;
#endif

const double EPSm3  = 1.0e-3;
const double EPSm4  = 1.0e-4;
const double EPSm5  = 1.0e-5;
const double EPSm6  = 1.0e-6;
const double EPSSP  = 6.0e-7;
const double EPSm7  = 1.0e-7;
const double EPSm8  = 1.0e-8;
const double EPSm12 = 1.0e-12;
const double EPSDP  = 6.0e-14;
const double EPSm14 = 1.0e-14;
const double EPSm20 = 1.0e-20;
const double EPSm30 = 1.0e-30;

const int StrMax    = STR_MAX;

enum lev {WARNING, ERROR, REMARK};

extern "C" {
  void       message (const char *routine, const char *txt, int level);

  double     dclock  ();
  float      sclock  ();

  double   *dvector  (integer nl, integer nh);
  double  **dmatrix  (integer rl, integer rh, integer cl, integer ch);
  double ***d3matrix (integer rl, integer rh, integer cl, integer ch,
		      integer dl, integer dh);

  float    *svector  (integer nl, integer nh);
  float   **smatrix  (integer rl, integer rh, integer cl, integer ch);
  float  ***s3matrix (integer rl, integer rh, integer cl, integer ch,
		      integer dl, integer dh);

  integer   *ivector  (integer nl, integer nh);
  integer  **imatrix  (integer rl, integer rh, integer cl, integer ch);
  integer ***i3matrix (integer rl, integer rh, integer cl, integer ch,
		       integer dl, integer dh);

  void freeDvector  (double    *v, integer nl);
  void freeDmatrix  (double   **m, integer nrl, integer ncl);
  void freeD3matrix (double  ***t, integer nrl, integer ncl, integer ndl);

  void freeSvector  (float     *v, integer nl);
  void freeSmatrix  (float    **m, integer nrl, integer ncl);
  void freeS3matrix (float   ***t, integer nrl, integer ncl, integer ndl);

  void freeIvector  (integer   *v, integer nl);
  void freeImatrix  (integer  **m, integer nrl, integer ncl);
  void freeI3matrix (integer ***t, integer nrl, integer ncl, integer ndl);
}

#endif

