##############################################################################
# Makefile for spectral element foundation routines, and shared definitions.
#
# $Id$
##############################################################################

ifdef DEBUG
  DEFINES += -DDEBUG
endif
ifdef MPI
  DEFINES += -DMPI
endif

CPPFLAGS = -I. -I$(SEM)/include
LIBDIR   = -L$(SEM)/lib/$(ARCH)

# ----------------------------------------------------------------------------
# Define architectures.
#
ARCH = $(shell uname -s)
MACH = $(shell uname -m)
ifeq ($(ARCH),Linux)
  ifeq ($(MACH),alpha)
    ARCH = Linux-alpha
  endif
  ifeq ($(MACH),i586)
    ARCH = Linux-iX86
  endif
  ifeq ($(MACH),i686)
    ARCH = Linux-iX86
  endif
  ifeq ($(MACH),x86_64)
    ARCH = Linux-x86_64
  endif
endif

# ----------------------------------------------------------------------------
# Compiler options for specific machines.
#

ifeq ($(ARCH),Linux-x86_64)             # -- Linux Opternon w GNU compilation.
CXX      = g++
LD       = g++
CC       = gcc
FC       = g77
ifdef DEBUG
  OPT    = -gstabs+ -Wno-deprecated
else
  OPT    = -O3 -w -Wno-deprecated
endif
CFLAGS   = $(DEFINES) $(OPT)
FFLAGS   = $(DEFINES) $(OPT)
CXXFLAGS = $(DEFINES) $(OPT)
LDFLAGS  = $(OPT) $(LIBDIR) -lfem -lvec $(F77LIBS) -lm
RANLIB   = ranlib
#F77LIBS  = -llapack -lblas -lg2c
F77LIBS  = -llapack -lgoto_opt64-r0.96 -lg2c
LIBDIR  += -L/usr/local/lib
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),Darwin)               # -- Mac OSX.
#DEFINES += -D_VECTOR_ARCH
CXX      = g++
LD       = g++
CC       = gcc
FC       = g77
AR	 = libtool
ARFLAGS  = -static -o
RANLIB   = ranlib
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3
endif
CFLAGS   = $(DEFINES) $(OPT)
FFLAGS   = $(DEFINES) $(OPT)
CXXFLAGS = $(DEFINES) $(OPT)
LIBDIR   += -L/sw/lib
LDFLAGS  = $(OPT) $(LIBDIR) -lfem -lvec $(F77LIBS)
F77LIBS  = -framework Accelerate -lg2c
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),SunOS)              	# -- Sun.
CXX    = CC
LD     = CC
CC     = cc
FC     = f77
ifdef DEBUG
  OPT    = -g -xarch=v9
else
  OPT    = -fast -xarch=v9
endif
CFLAGS   = $(DEFINES) $(OPT)
FFLAGS   = $(DEFINES) $(OPT)
CXXFLAGS = $(DEFINES) $(OPT) -library=iostream
LDFLAGS  = $(LIBDIR)  $(OPT) -library=iostream
CLDFLAGS = $(LIBDIR)  $(OPT)
RANLIB   = 
F77LIBS  = -xlic_lib=sunperf
ifdef MPI
  MPIHOME   = /opt/lam-6.3.3b20
  CPPFLAGS += -I$(MPIHOME)/include
  LIBDIR   += -L$(MPIHOME)/lib
  LDFLAGS  += -lfem_mp -lvec $(F77LIBS) -lmpi -llam -lsocket
else
  LDFLAGS  += -lfem -lvec $(F77LIBS)
  CLDFLAGS += -lfem -lvec $(F77LIBS)
endif
endif

ifeq ($(ARCH),Linux-alpha)              # -- Linux alpha w Compaq C & f90.
CXX    = g++			        # Use g++ as C++ and linker.
LD     = g++			        # This avoids I/O related problems.
CC     = ccc
FC     = fort
ifdef DEBUG
  GOPT    = -gstabs+
  COPT    = -g2
else
  GOPT    = -O4
  COPT    = -fast -tune host
endif
CFLAGS   = $(DEFINES) $(COPT)
FFLAGS   = $(DEFINES) $(COPT)
CXXFLAGS = $(DEFINES) $(GOPT)
LDFLAGS  = $(LIBDIR)  $(GOPT)
RANLIB   = ranlib
LIBDIR  += -L/usr/local/lib
F77LIBS  = -lcxml -lfor -lots
ifdef STATIC
  LDFLAGS += -non_shared
endif
ifdef NOCXML     # -- For testing, avoid Compaq math libs.
  F77LIBS =  -llapack -lgemm -lblas -lots -lfor -lg2c
endif
ifdef MPI
  MPIHOME   = /usr/local/mpich
  LDFLAGS  += -lfem_mp -lvec $(F77LIBS) -lcpml \
              -lmpich -lsnet_vipl -lvipl -ldl
# MPIHOME   = /usr/local/lam
# LDFLAGS  += -lfem_mp -lvec $(F77LIBS) -lcpml \
#             -lmpi -ltrillium -ltstdio -lt -largs
  CPPFLAGS += -I$(MPIHOME)/include
  LIBDIR   += -L$(MPIHOME)/lib
else
  LDFLAGS += -lfem -lvec $(F77LIBS) -lcpml
endif
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),Linux-iX86)               # -- Linux Intel w GNU compilation.
CXX      = g++
LD       = g++
CC       = gcc
FC       = g77
ifdef DEBUG
  OPT    = -gstabs+ -Wno-deprecated
else
  OPT    = -O3 -w -Wno-deprecated
endif
CFLAGS   = $(DEFINES) $(OPT)
FFLAGS   = $(DEFINES) $(OPT)
CXXFLAGS = $(DEFINES) $(OPT)
LDFLAGS  = $(OPT) $(LIBDIR) -lfem -lvec $(F77LIBS) -lm
RANLIB   = ranlib
F77LIBS  = -llapack -lblas -lg2c
#F77LIBS  = -llapack -lfastblas -lg2c
LIBDIR  += -L/usr/local/lib
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),IRIX)			# -- SGI 32-bit systems.
CXX	 = CC
LD	 = C
RANLIB   = 
CFLAGS   = $(DEFINES) -fullwarn -O2
FFLAGS   = $(DEFINES) -O2
CXXFLAGS = $(DEFINES) -fullwarn -O2
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = $(LIBDIR) -lfem -lvec $(F77LIBS) -lmalloc -lm
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),IRIX64)			# -- SGI 64-bit systems.
CXX	 = CC
LD	 = CC
CC       = cc
FC       = f77
CFLAGS   = $(DEFINES) -n32 -mips4 $(OPT)
FFLAGS   = $(DEFINES) -n32 -mips4 $(OPT)
CXXFLAGS = $(DEFINES) -n32 -mips4 $(OPT)
LDFLAGS  = -n32 -mips4 $(OPT) $(LIBDIR) -lfem -lvec $(F77LIBS) -lm
RANLIB   =
ifdef DEBUG
  OPT    = -g
else
  OPT    =  -O3
endif
F77LIBS  = -lcomplib.sgimath -lftn
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),OSF1)			# -- DEC OSF1 Alpha.
CXX	 = cxx
LD	 = cxx
FC       = f90
CC       = cc
CFLAGS   = $(DEFINES) -fast -tune host
FFLAGS   = $(DEFINES) -fast -tune host
#CXXFLAGS = $(DEFINES) -fast -tune host -strict-ansi -D__USE_STD_IOSTREAM
CXXFLAGS = $(DEFINES) -fast -tune host
RANLIB   =
# -- The EV5 versions of the BLAS, LAPACK are quicker than the EV6's!
FORTLIBS = /opt/EV5/usr/opt/XMDLIB5500/libcxml_ev5.a -lfor
#FORTLIBS = -ldxml -lfor
# -- This redefinition is required on Compaq machines:
#YACC     = bison -o y.tab.c
ifdef MPI
  MPIHOME   = /usr/local/mpi
  CPPFLAGS += -I$(MPIHOME)/include  
  LIBDIR   += -L$(MPIHOME)/lib/alpha/ch_shmem
  LDFLAGS   = $(LIBDIR) -lfem_mp -lvec $(FORTLIBS) -lm -lmpi
# -- This is needed for Compaq MPI:
# LDFLAGS  += -lelan
else
  LDFLAGS   = $(LIBDIR) -lfem -lvec $(FORTLIBS) -lm
endif
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),UNIX_System_V)    	# -- Fujitsu VPP300.
DEFINES += -D_VECTOR_ARCH
CC       = vcc
FC       = frt
CXX      = /home/home01/565/hmb565/C++/bin/CC
LD       = $(CXX)
OPT      = -Wv,-m3,-Om,-te -Kpopt
CFLAGS   = $(DEFINES) -Wv,-m3,-Om -Kpopt
FFLAGS   = $(DEFINES) -Wv,-m3,-te
CXXFLAGS = $(DEFINES)  $(OPT)
RANLIB   =
ifdef MPI
  CPPFLAGS += -I/usr/lang/mpi/include
  LIBDIR   += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
  LIBDIR   += -L/usr/lang/mpi/lib
  LIBS      = -lfem_mp -lvec
  LIBS     += -llapack -lblasvpp -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
  LIBS     += -lmpi -lmp
  OPT      += -Wl,-P
else
  LIBDIR   += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
  LIBS      = -lfem -lvec
  LIBS     += -llapack -lblasvpp -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
endif
LDFLAGS  = $(OPT) $(CPPFLAGS) $(LIBDIR) $(LIBS)
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),SUPER-UX)         	# -- NEC SX-4/5.
DEFINES += -D_BUILTIN_ -D_SX -D_VECTOR_ARCH
CXX       = c++
LD        = c++
FC        = f90
CFLAGS    = $(DEFINES) -pvctl,fullmsg,loopcnt=20000 -Onooverlap \
           -hfloat0 -hnostkchk
FFLAGS    = $(DEFINES) -Wf'-pvctl,fullmsg -ptr byte -b' -float0 -ew
CXXFLAGS  = -V $(DEFINES) -pvctl,loopcnt=20000,fullmsg -acct $(OPT)
RANLIB    =
CPPFLAGS += -I/usr/include/CC -I/usr/local/include/CC
LIBDIR   += -L/usr/local/lib
OPT       = -O nooverlap -Nstkchk
NOVEC     = -c -Nstkchk -acct -dir novec
FORTLIBS  = -llapack_64 -lblas_64 -lfft_64
ifdef MPI
  LIBS    = -lfem_mp -lvec $(FORTLIBS) -lm -lC -lmpi -lpthread
else
  LIBS    = -lfem -lvec $(FORTLIBS) -lm -lC
endif
LDFLAGS   =  -verbose -V -f90libew $(CXXFLAGS) $(CPPFLAGS) $(LIBDIR) $(LIBS)
# -- Add this to enable profiling:
#LDFLAGS += -p
CLDFLAGS  = $(LDFLAGS)
endif

# ----------------------------------------------------------------------------
# Installation of header files.
#
default:
	cp -f sem.h analysis.h auxfield.h bcmgr.h boundary.h bsys.h	\
	condition.h data2df.h domain.h edge.h element.h family.h feml.h	\
	field.h flowrate.h geometry.h history.h integration.h matrix.h 	\
	mesh.h misc.h numbersys.h particle.h pressure.h statistics.h 	\
	../include

#-----------------------------------------------------------------------------
# Force installation of headers from another make.
#
sem:
	cd $(VPATH); $(MAKE)

# -----------------------------------------------------------------------------
# Clean up.
#
clean:
	rm -f *.o *.O *.a *~ *.ti SunWS_cache
	rm -rf ILDUMPS ii_files ptrepository
