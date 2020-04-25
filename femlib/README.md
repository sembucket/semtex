# Femlib README {#femlib_readme}

Femlib README
=============

Synopsis
--------

Routines to support finite/spectral element operations, Fourier
transforms, and memory exchanges.  The routines are grouped together
to create the static object library libfem.a.

Notes
-----

1. All inter-process memory exchange operations in semtex are handled
in file message.c.  Memory exchanges are handled by MPI calls, which
are quarantined inside message.c.  For single-process compilation, all
message-passing routines in message.c become empty stubs.

2. Lowest-level access to orthogonal polynomial routines is provided
in files polyops.c and polylib.c, which have a long history dating
back to Einar Ronquist's equivalent Fortran routines of the early
1990s, followed by Ron Henderson's re-implementation in C, and
subsequent tidying-up, extensions and bug-fixing by Spencer Sherwin
(in file polylib.c, which is adopted from Nektar/Nektar++.  In fact
only a few routines from polylib.c are used, where the available
functionality was not available/correct in polyops.c.  Probably we
should consolidate/standardize to use only the 1D orthogonal
polynomial routines of polylib.c).

3. Various higher-level routines which provide access to the
polynomial routines of polyops.c and polylib.c are supplied in
operators.c, the latter largely directed at providing 2D operators (as
opposed to 1D operators of polylib.c and polyops.c).

4. Because semtex is essentially a 2D spectral element x Fourier code,
discrete Fourier transform routines are also supplied by femlib, in
file fourier.c, which provides routines to compute multiple 1D
real-complex FFT/IFFT.  All the lower-level implementation for
discrete Fourier transforms is in Fortran77.  The default FFT is now a
2,3,5 prime factor complex-complex FFT written by Clive Temperton;
this is provided in file temfft.F, along with top-level pack/unpack
interface routine which converts the complex-complex FFT to a
real-complex FFT.  For debugging purposes, there are also (somewhat
slower) real-complex FFT routines from netlib/FFTPACK, provided in
netlib.f; to use these routines, compile fourier.c with DEBUG_FFT
defined.  (An advantage of the FFTPACK routines is that they can be
used for arbitrary length FFTs.)  Within fourier.c there is also the
possibility to use a NEC-supplied FFT library routines if these are
available (_SX must be defined); this interface amounts to legacy code
but might be useful if you have an NEC machine and its (fast!) vendor
libraries.  Finally, radix-2 FFT routines from the back of Canuto et
al. (1988) are available too, in file canfft.f.

5. A short routine which supplies mapping data from standard 2D
elemental layout to 1D vectors is given in mapping.c.

6. Family.c provides routines to store matching 1D storage (somewhat
like the operation of smart pointers).  Largely replaced by
src/family.cpp and smart pointers.  Needs to be reviewed - could C++
smart pointers now do all of this?

7. Filter.c.  This provides a routine to compute coefficients of
smooth lowpass filter (the "Boyd-Vandeven" filter).

8. Femlib provides a yacc-based function parser, which is much used
within semtex to maintain a set of tokens and to parse functions.  The
code for this was originally based on hoc3 in "The UNIX programming
environment" by Kernighan and Pike, and is contained in file
initial.y.

9. Header files.

* The principal header file to allow use of femlib library routines
  from C++ is femlib.h, which makes use of static member functions to
  provide name-qualified access to the routines in the various C and
  F77 files that comprise femlib.  For example Femlib::value() allows
  a clearly identifiable interface to routine yy_interpret() which is
  part of the built-in function parser.  (This naming mechanism was
  incorporated into semtex before namespaces were adopted in C++.
  Nowadays I'd probably use namespaces instead.)

* If you want to call femlib routines from C instead of C++ (as is
  done for example in various utility routines), you can use the
  prototypes provided in cfemlib.h.

* The standard semtex real and integer data types real_t and int_t are
  given by typedefs in cfemdef.h, which supplies the lowest-level
  definitions and macros used in semtex.

* Various pre-defined TOKEN values are supplied to the function
  parser.  These are contained in file defaults.h.

______________________________________________________________________________
$Id: README.md,v 9.1 2020/01/06 04:35:44 hmb Exp $
