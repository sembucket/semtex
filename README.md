# Semtex top-level README file (semtex/README)

# List of other sub-pages

* Page [semtex/veclib/README]  (@ref veclib_readme)
* Page [semtex/femlib/README]  (@ref femlib_readme)
* Page [semtex/utility/README] (@ref utility_readme)

_____________________________________________________________________________

semtex README
=============

Synopsis
--------

The semtex package is for spectral element solution of (incompressible
Navier--Stokes) time-varying and elliptic PDEs.  Geometries
accommodated are 2D Cartesian or cylindrical, and/or 3D
(a.k.a. "2+1/2D", i.e. 2D X Fourier), again Cartesian or cylindrical.
In fact, the main emphasis is on direct numerical simulation of
incompressible flows, with the elliptic solver (which is used by the
DNS solver) available as a stand-alone adjunct.

Flows can have 2 or 3 velocity components, so in terms of number of
spatial dimensions (D) and velocity components (C), the available
combinations are: 2D2C, 2D3C, 3D3C. If desired, a scalar can be added
to DNS (and also, if requested, advected in a frozen/supplied velocity
field).

Spectral elements are used to discretise planar geometries, with
solution variations in the third (homogeneous/periodic) direction
accommodated via Fourier expansions (i.e. spatial discretisation is 2D
spectral element X Fourier).

References
----------

The first of these is the recommended starting point; it provides an
introductory overview of the code and its utilities.  The second gives
details of the cylindrical-coordinate formulation.

1. Blackburn HM, Lee D, Albrecht T & Singh, J (2019) Semtex: a
spectral element–Fourier solver for the incompressible Navier–Stokes
equations in cylindrical or Cartesian coordinates. Computer Physics
Communications V245: 106804-1–13.

2. Blackburn HM & Sherwin SJ (2004) Formulation of a Galerkin spectral
element–Fourier method for three-dimensional incompressible flows in
cylindrical geometries. Journal of Computational Physics V197N2:
759–778.

Technical details
-----------------

The code is designed to compile and run on Unix-based operating
systems including Linux and OSX (and has been successfully compiled
and run on a wide variety of such systems for over two decades).
Optional multi-process parallelisation across Fourier modes is
supported via MPI.

The method uses continuous-Galerkin Legendre--Gauss--Lobatto
nodal-based spectral-element elliptic equation solvers, a 2,3,5
prime-factor FFT, and a backward-time (a.k.a. stiffly-stable)
equal-order (P_N-P_N) velocity-pressure semi-explicit timesplitting
method for DNS, with first-, second-, or third-order time integation.
Skew-symmetric forms of the nonlinear terms in the NSE are the default
for robust operation. No form of dealiasing is used when creating
nonlinear product terms.

Elliptic equations are by default solved using direct methods
(Cholesky decomposition allied with element-level static
condensation), but optionally can be solved using (Jacobi)
preconditioned conjugate-gradient methods allied with tensor-product
methods.  The code typically operates most efficiently with moderate
spectral element polynomial orders (4-10).

The standard real data type in semtex is 64-bit double precision
(double in C, DOUBLE PRECISION in Fortran) while the standard integer
type is usually 32-bit (typically, whatever int/INTEGER is on the
machine you are running on).  These two types are referred to by the
semtex-standard typedefs real_t and int_t (see femlib/cfemdef.h).
Strings are by default 2048 (STR_MAX) characters long.

Distribution
------------

The src directory contains source for the central classes of spectral
element solvers.  Two C/F77 libraries called veclib and femlib (source
files in eponymous directories) contain low-level algebraic,
polynomial, string parsing and message-passing routines.  The
upper-level codes are written in C++.

Source for application programs can be found in three upper-level
directories:

elliptic:   elliptic (Laplace, Poisson, Helmholtz) solver.  
dns:        Navier--Stokes DNS solver (uses same elliptic solver routines).  
utility/*:  various utility programs.

A user guide and HTML documentation is provided in the doc directory.

Required third-party software
-----------------------------

1. A Unix-based operating system (which mostly can readily supply the
   following things, either as-is or via package installations):

2. make (Gnu's version, which is usually now the standard supplied code);

3. C++, C and Fortran compilers (OSX needs both vendor-supplied Xcode
   for C/C++ and a 3rd-party Fortran compiler, installed e.g. via
   homebrew/macports/fink as part of the gcc compiler suite);

4. LAPACK and BLAS libraries;

5. bison or yacc.

Optional/useful third-party software
------------------------------------

1. cmake: baseline build method, but make can alternatively be used; see below.

2. (Open) MPI and associated header files: to make optional parallel solvers.

3. Tecplot and/or Paraview, VisIt: for visualisation and some postprocessing.

4. Supermongo: an old 2D plotting package, for which some macros are
supplied.  You may need to buy a licence for supermongo.

Preliminary notes for Mac OS X users
------------------------------------

The code is designed to compile on Unix machines, including Mac OS X
(which is based on BSD Unix) and of course, Linux. In fact, since
2004, Semtex has predominantly been developed on OS X.  For OS X, you
will also need to have installed: Xcode (from Apple), and (at least) a
Fortran 77 (or F90, F95) compiler.  For the latter, it is usual to
install one of the GNU/Unix gcc compiler suites, available through one
of the standard open-source software ports for OS X (macports,
homebrew, or fink).  You might also wish to install an MPI setup, such
as openmpi.  There is no need to install your own BLAS or LAPACK, as
these usually come as a standard part of Xcode (in the Accelerate
package).

Building - Introduction
-----------------------

Semtex presently supports two alternative building methods:

1. Out-of-source compilation using cmake (this is the recommended
approach since it is most likely to succeed out-of-the-box: failure
is usually associated with non-existence of some required 3rd-party
software, which will be reported by cmake);

2. In-source compilation using make (this is a little more complicated
but more readily allows fine-scale control of compilers, choice of
BLAS routines, and the like).  If your machine is old and does not
have cmake2.8 or above, in-source building may be your only option.

Note that the use of vendor-supplied BLAS routines such those in
Intel's MKL or AMD's ACML may quite substantially improve the
performance of the DNS code.  Most of the improvement comes from
having a well-optimised version of DGEMM (matrix-matrix multiply),
which is heavily used by dns.  Search GotoBLAS or OpenBLAS.  On OSX,
you get this performance by default via the Xcode/Accelerate package;
on Linux, try using ACML or MKL, depending on machine architecture.
Speed improvements to be had by choosing a good BLAS version usually
far outweigh what can be achieved via choice of compiler or compiler
optimisation flags.

Out-of-source building with cmake
---------------------------------

In this case, all the supplied directives are in files called
CMakeLists.txt in various directories (and you can ignore any supplied
Makefiles).

Make a build directory at the top level of the directory tree (which
contains this README file) as build_dir: the name used to replace
build_dir is your choice (we usually use "build").  This is where
all the executables will end up, and you could conveniently add this
directory to your PATH.  If you do not have MPI/openmpi installed on
your system, only the serial versions of the solvers dns and elliptic
will be built.  If cmake can find a working MPI installation, parallel
versions will also be built (which are called dns_mp and elliptic_mp).
All the other resulting files are for serial execution only.

  %> mkdir build_dir
  %> cd build_dir
  %> cmake ..
  %> make
  %> ctest

You should find that all test regression tests are reported as passed.
Because these tests exercise the solvers as well as a range of
utilities, SUCCESSFUL COMPLETION OF THE CTEST STAGE INDICATES THAT YOU
HAVE A WORKING SET OF EXECUTABLES.

You can also explicitly turn off MPI-based parts of compilation,
and/or make executable versions for debugging if required, via
command-line flags to cmake, e.g.: cmake -DWITH_MPI=OFF -DDEBUG=ON ..

In-source building with make
----------------------------

In this case, the supplied Makefiles are relevant (and CMakeLists.txt
files can be ignored).  System-dependent compilation flags can be
found/edited in src/Makefile -- you will likely find a set of flags
that either suit your machine or can be edited as required. The
executables will be scattered around in their various source
directories, in which case you may soft-link them into some central
directory (e.g. ~/bin) to get them in your PATH.

From the top-level directory (contains this file, README), do

  %> make test

That will build libraries veclib.a, femlib.a and a few core
executables (enumerate, compare, and dns), then run regression tests
on dns.  If everything is fine, compilation will proceed without a
hitch and at the end, various tests will run and report as passed.
Again, successful conclusion of this stage that you have a working set
of (in this case, serial-only) executables.

If compilation does not complete then the likely alternatives are that
(a) you are missing some required 3rd-party software or (b)
compiling/linking options are incorrect for your machine.

If it is (a) and you are on a computer that uses modules (try "module
list"), you may be able to cure the problem by loading a correct set
of modules (try "module avail").  Otherwise you may have to install
BLAS/LAPACK libraries with whatever package tool exists on the
machine, or in the worst case by compiling and installing these for
yourself.  If you are using OSX, you will need to install Xcode, the
Xcode command-line tools, and install a Fortran compiler either as a
stand-alone package or by installing the gnu compiler suite using one
of fink, homebrew or macports.

If it is (b), the compiler and link options are all set in
src/Makefile; after running "uname -s; uname -m; uname -n" you should
be able to find the set of flags/link options being used in
src/Makefile.  Often it will be "fairly" obvious which flags are in
error -- just how obvious probably depends on your experience with the
machine in question (try some googling if you have no clue).  Try
compiling both the veclib and femlib libraries by running make in
their directories (this will also show you the compiler flags being
used).  If the libraries can be compiled and archived (to produce
libvec.a and libfem.a) the most likely issue is the final linking of
the executables for dns, enumerate and compare.  Once you can get the
libraries to compile "by hand", go back up to this (README) directory
and do "make clean; make libs".  When that completes OK, go to the dns
directory and do "make clean; make", with the aim of compiling and
linking the executable "dns".  If dns is not produced, you will need
to change the linking options in src/Makefile.  Check and edit the
link flags chosen from src/Makefile and those which work on your
machine until you are able to compile all of libfem.a, libvec.a and
dns.  At this point, go back and do "make test" as above - very likely
the two required serial utilities will also now compile.  When "make
test" completes without error, try compiling all the remaining serial
utilities:

  >% make all

At this stage, and if you have MPI, you may wish to try building and
linking the parallel library libfem_mp.a and solver dns_mp (typically
run like this: "mpirun -np <n_cores> dns_mp").

  >% make parallel

You should find executables "dns" and "dns_mp" now appear in the dns
directory.  Congratulations - you are finished!  Putting all the
executables in your PATH could be your next task before going on to do
some computing.  Given all the above description, you could be
forgiven for trying the cmake route first.

Files types used by semtex
--------------------------

Semtex uses a base input file which describes the mesh and boundary
conditions.  We call this a "session" file and typically it has no
root extension.  It is written in a format patterned on HTML, which we
have called FEML (for Finite Element Markup Language).  There are a
number of example session files in the mesh directory.  Other files
have standard extensions:

session.num  Global node numbers, produced by enumerate utility.  
session.fld  Solution/field file.  Binary format by default.  
session.rst  Restart file.  Read in to initialise solution if present.  
session.avg  Average file.  Used to store time-mean averages.

After writing a new session file it is best to run meshpr on it before
trying to use it for simulations.  Meshpr will catch most of the
easier-to-make errors.  You can also plot up the results using
SuperMongo or other utility as a visual check.  The rectmesh utility
can be used to produce session files for simple rectangular domains -
you can edit their tokens and boundary conditions as you see fit.

Utilities
---------

Can be found in the utility directory, including:

enumerate Generate global node numbering, with RCM optimization.  
compare   Generate restart files, compare solutions to a function.  
convert   Convert field file formats (IEEE-big/little, ASCII).  
meshpr    Generate 2D mesh locations for plotting or checking.  
addfield  Add vorticity, also divergence, helicity, etc to a field file.  
sem2tec   Convert field files to Tecplot format (OK with Paraview/VisIt too).  
project   Convert a field file to a different order interpolation.  
interp    Interpolate field file from one 2D mesh to another.  
probe     Probe field file at specified points.  
noiz      Add a random perturbation to a field file.  
calc      An interactive calculator that calls the built-in function parser.  
rstress   Compute a Reynolds stresses from a .avg file, subtract averages.  
rectmesh  Generate a start-out session file from a list of x and y values.  

Most executables have a -h command line option which gives a help
prompt.  If that is insufficient help, please read the header section
of the associated source files or the semtex user guide.

User guide and doxygen documentation
------------------------------------

Assuming you have a latex/pdflatex system installed: 

  %> cd doc; make

to produce userguide.pdf.

Assuming you have doxygen (and dot) also installed:

  %> cd doc; make doxygen

and open doc/html/index.html in a browser.

Author
------

Hugh M. Blackburn  
Department of Mechanical & Aerospace Engineering  
Monash University  
Vic 3800  
Australia  
mailto:hugh.blackburn@monash.edu

Licencing
---------

Please note that semtex is a research code and is provided 'as-is'
with the understanding that it will be used mainly by other
computer-literate researchers in computational fluid dynamics.  It is
not guaranteed to work or to provide correct results, and neither the
author or any employer accept any liability for detriment or loss
consequent on your use of the code.  Source code for semtex is free
software as released for use under the GNU GPL2 public license; please
consult the accompanying COPYING file for details.

List of major revisions
-----------------------

1995: Semtex-1:   ase 2D version of code completed in C.         
1996: Semtex-2:   C --> C++ conversion completed.  
1996: Semtex-3:   2D Cartesian/Fourier (i.e. 3D periodic) spaces.  
1997: (Scat, separate code) Heat/scalar transfer supported.  
1997: Semtex-4.0: Cylindrical solutions in 2D or 3D supported.      
1997: Semtex-4.1: Generalized prime factor FFT routines.  
1997: Semtex-4.2: Massless particle tracking.  
1997: Semtex-5.0: Concurrent execution with MPI.  
1998: Semtex-5.1: Improved vectorization & IO performance.  
1999: Semtex-5.3: Mixed/Robin BC type added.  
2003: Semtex-5.5: Adopt standard C++ libraries wherever possible.       
2004: Semtex-6:   Cylindrical coordinate/3D code exponentially convergent.  
2004: Semtex-6:   Mac OSX port -- filenames no longer case-sensitive.  
2010: Semtex-7:  "I do not now recall".  
2016: Semtex-8:   Generalised body forces supported in DNS.  
2018: Semtex-9:   Scalar transport and DNS code merged.  
2019: Semtex-9.1: cmake adopted as baseline compilation system.          

------------------------------------------------------------------------------

$Id: README.md,v 1.1 2020/01/06 04:35:44 hmb Exp $
