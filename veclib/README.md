# Veclib README {#veclib_readme}

Veclib README
=============

Synopsis
--------

Vector algebra primitives.  The routines are linked together into the
static library libvec.a.

Notes 
-----

1. The library descends from the Intel iPSC/2 VECLIB library (circa
1990).  Available vector operations complement the BLAS.  See the
iPSC/2 Programmer's Reference Manual and the C VECLIB Routine Summary
for further information.

2. Not all the original VECLIB routines are implemented: for example
those which are already in the BLAS are omitted and there are some
additions, notably those for memory management.  (The BLAS are given
C-compatible call interfaces by macros in cveclib.h.)

3. The names of the associated C files in this directory usually start
with "x" indicating they implement the related routines in a variety
of precisions (see list below).

4. The naming convention for triad and some of the miscellaneous
routines follows an abbreviated RPN-like syntax based on the letters
"s" (scalar), "v" (vector), "p" (plus), "m" (minus), "t" (times) .
For example, "svmvt" stands for "scalar,vector,minus,vector,times",
while "svvtp" stands for "scalar,vector,vector,times,plus".  Check the
corresponding files (xsvmvt.c and xsvvtp.c) to see how the strings
correspond to operations.

5. No checks are made for length of vectors, negative lengths.
	
6. Negative increments are allowed; behaviour is the same as for the
BLAS, i.e. if the vector increment is negative, the supplied
vector address is the last item to be accessed, with vector traverse
beginning at the supplied address + (n-1) * abs (inc), (further along
the vector).

7. Header/function prototype files:

* Static member functions in veclib.h make the various routines herein
  available for C++ routines with the preceding name Veclib, along
  with use of overloading e.g. Veclib::copy() might resolve to
  icopy(), scopy() or dcopy from file xcopy.c, depending on the type
  of arguments supplied.  This use of static member functions would
  nowadays be replaced using C++ namespaces (which were not part of
  standard C++ when semtex was first written).  Of course, one could
  make the approriate calls directly to the underlying C routines, but
  this way it is always clear from which library the routines derive.
  The corresponding prototypes to allow use of veclib routines from C
  are supplied in cveclib.h.

* Equivalent static member functions defined in blas.h and lapack.h
  allow "namespace" functionality with BLAS and LAPACK calls (along
  with overloading).  For example Blas::dot() might resolve to sdot_()
  or ddot_().  These header files aren't really part of veclib, but
  are supplied here for want of an alternative location.  (The
  trailing underscore on these names is typical of how routines
  compiled from Fortran source will appear in UNIX object files: the
  macro F77NAME from femlib/cfemdef.h automates the renaming.)

* The file utility.h provides low-level definitions and templates as
  well as prototypes for small routines in util.c.  One low-level
  routine from here that is often used throughout the code base
  (**without** the preceding Veclib::!) is message().

<pre>

Key:
----
#:	not implemented.
%:	replaced by BLAS routine.
x:	implemented with indicated precision.
*:	an extension to VECLIB.
-:	not applicable.

-----------------------------------------------------------------------------
Section         Name               Available Precision	d   s   c   z   i
-----------------------------------------------------------------------------

UTILITIES:
*               message
*               efopen
*               printxvector                            x   x           x

MEMORY MANAGEMENT:
*               vector  / free_vector.................. x   x   x   x   x
*               matrix  / free_matrix.................. x   x   x   x   x
*               3matrix / free_3matrix................. x   x   x   x   x

MATHEMATICAL PRIMITIVES:
%               swap
                copy................................... x   x           x
                fill................................... x   x           x
                neg.................................... x   x           x
                vneg................................... x   x           x
*               znan................................... x   x
*               vsgn................................... x   x           x
                sadd................................... x   x           x
                vadd................................... x   x           x
                vsub................................... x   x           x
                smul................................... x   x           x
                vmul................................... x   x           x
%               scal
                sdiv................................... x   x           x
                vrecp.................................. x   x
                vdiv................................... x   x           x
*               zero................................... x   x           x
*               spow................................... x   x

OTHER MATHEMATICAL FUNCTIONS:
                vabs................................... x   x           x
#               vmax
#               vmin
                vamax.................................. x   x           x
#               vamin
                vpow................................... x   x
                vexp................................... x   x
                vlg10.................................. x   x
                vlog................................... x   x
                vatan.................................. x   x
                vatn2.................................. x   x
                vcos................................... x   x
                vsin................................... x   x
                vsqrt.................................. x   x
                random (ranu).......................... x   x
                vrandom................................ x   x
*               normal................................. x   x
*               vnormal................................ x   x
*               fft.................................... x   x   x   x
*               vhypot................................. x   x
*               vtanh.................................. x   x

TRIAD OPERATIONS:
%               axpy
                svmvt.................................. x   x
                svpvt.................................. x   x
                svtsp.................................. x   x
                svtvm.................................. x   x
                svtvp.................................. x   x
                svvmt.................................. x   x
                svvpt.................................. x   x
                svvtm.................................. x   x
                svvtp.................................. x   x
*               svvtt.................................. x   x
                vvmvt.................................. x   x
                vvpvt.................................. x   x
                vvtvm.................................. x   x
                vvtvp.................................. x   x
                vvvtm.................................. x   x
*               vvvtt.................................. x   x

RELATIONAL PRIMITIVES:
#               eq
                seq....................................                 x
#               ge
                sge.................................... x   x           x
#               gt
#               sgt
                sle.................................... x   x           x
                slt.................................... x   x           x
#               ne
                sne.................................... x   x           x

LOGICAL PRIMITIVES:
#               land
#               lnot
#               lor
#               lsand
#               lsor
#               lcopy
#               lfill

REDUCTION FUNCTIONS:
%               asum
%               dzasum
%               scasum
                sum.................................... x   x           x
%               nrm2
%               dot
%               dotc
%               dotu
                i_max.................................. x   x           x
                i_min.................................. x   x           x
                icount................................. -   -   -   -   x
                ifirst................................. -   -   -   -   x
#               ilast
                lany................................... -   -   -   -   x
*               lxsame................................. x   x           x

CONVERSION PRIMITIVES:
#               vcmplx
#               vconjg
                vdble..................................     x
                vsngl.................................. x
                vfloa.................................. x   x   -   -   -  
*               brev................................... x   x           x

MISCELLANEOUS FUNCTIONS:
                scatr.................................. x   x           x
                gathr.................................. x   x           x
*               gathr_scatr............................ x   x           x
*               scatr_sum.............................. x   x           x
*               gathr_sum.............................. x   x           x
*               gathr_scatr_sum........................ x   x           x
                ramp................................... x   x           x
                clip................................... x   x           x
*               clipup................................. x   x           x
*               clipdn................................. x   x           x
                iclip.................................. x   x           x
                cndst.................................. x   x           x
                mask................................... x   x           x
#               fft
#               ifft
%               rot
%               rotg
#               folr
#               solr
#               lbidi
#               trfac
#               ubidi
                vpoly.................................. x   x
*               polint................................. x   x
*               splint................................. x   x
*               mxm.................................... x   x
*               mxv.................................... x   x
*               mxva................................... x   x
*               vvtvvtp................................ x   xw
*               vvtvvtm................................ x   xw
*               svvttvp................................ x   x

</pre>

$Id: README.md,v 1.1 2020/01/06 04:35:45 hmb Exp $
