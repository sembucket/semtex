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
endif

# ----------------------------------------------------------------------------
# Compiler options for specific machines.
#

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
  LDFLAGS  += -lfem_mp -lalp $(F77LIBS) -lmpi -llam -lsocket
else
  LDFLAGS  += -lfem -lalp $(F77LIBS)
  CLDFLAGS += -lfem -lalp $(F77LIBS)
endif
endif

ifeq ($(ARCH),Linux-alpha)              # -- Linux alpha w Compaq C & f90.
DEFINES += -D_Linux_alpha
CXX    = g++
LD     = g++
CC     = ccc
FC     = fort
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3
endif
CFLAGS   = $(DEFINES) -arch host -tune host -fast -O4
FFLAGS   = $(DEFINES) -arch host -tune host -fast -O4
CXXFLAGS = $(DEFINES) $(OPT)
LDFLAGS  = $(LIBDIR)
RANLIB   = ranlib
LIBDIR  += -L/usr/local/lib
ifdef STATIC
  F77LIBS  = -lcxml -lfor -lFutil -lUfor -lots -lg2c
  LDFLAGS += -static
else
  F77LIBS  = -lcxml -lfor -lots -lg2c
endif
ifdef NOCXML     # -- For testing, avoid Compaq math libs.
  F77LIBS =  -llapack -lgemm -lblas -lots -lfor -lg2c  
endif
ifdef MPI
  MPIHOME   = /usr/local/mpich
  LDFLAGS  += -lfem_mp -lalp $(F77LIBS) -lcpml \
              -lmpich -lsnet_vipl -lvipl -ldl
# MPIHOME   = /usr/local/lam
# LDFLAGS  += -lfem_mp -lalp $(F77LIBS) -lcpml \
#             -lmpi -ltrillium -ltstdio -lt -largs
  CPPFLAGS += -I$(MPIHOME)/include
  LIBDIR   += -L$(MPIHOME)/lib
else
  LDFLAGS += -lfem -lalp $(F77LIBS) -lcpml
endif
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),Linux-iX86)               # -- Linux Intel w GNU compilation.
CXX      = g++
LD       = g++
CC       = gcc
FC       = g77
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3 -w
endif
CFLAGS   = $(DEFINES) $(OPT)
FFLAGS   = $(DEFINES) $(OPT)
CXXFLAGS = $(DEFINES) $(OPT)
LDFLAGS  = $(OPT) $(LIBDIR) -lfem -lalp $(F77LIBS) -lm
RANLIB   = ranlib
F77LIBS  = -llapack -lblas -lg2c
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
LDFLAGS  = $(LIBDIR) -lfem -lalp $(F77LIBS) -lmalloc -lm
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
LDFLAGS  = -n32 -mips4 $(OPT) $(LIBDIR) -lfem -lalp $(F77LIBS) -lm
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
CFLAGS   = $(DEFINES) -fast -inline speed -fp_reorder \
           -speculate all -tune host -O4
FFLAGS   = $(DEFINES) -fast -tune host -O5
CXXFLAGS = $(DEFINES) -tune host -O2
RANLIB   =
FORTLIBS = -ldxml -lfor
ifdef MPI
  MPIHOME   = /usr/local/mpi/
  CPPFLAGS += -I$(MPIHOME)/include  
  LIBDIR   += -L$(MPIHOME)/lib/alpha/ch_shmem
  LDFLAGS   = $(CXXFLAGS) $(LIBDIR) -lfem_mp -lalp $(FORTLIBS) -lm -lmpi
else
  LDFLAGS   = $(CXXFLAGS) $(LIBDIR) -lfem -lalp $(FORTLIBS) -lm
endif
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),UNIX_System_V)    	# -- Fujitsu VPP300.
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
  LIBS      = -lfem_mp -lalp
  LIBS     += -llapack -lblasvpp -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
  LIBS     += -lmpi -lmp
  OPT      += -Wl,-P
else
  LIBDIR   += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
  LIBS      = -lfem -lalp
  LIBS     += -llapack -lblasvpp -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
endif
LDFLAGS  = $(OPT) $(CPPFLAGS) $(LIBDIR) $(LIBS)
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),SUPER-UX)         	# -- NEC SX-4.
DEFINES += -D_BUILTIN_ -D_SX
CXX       = c++
LD        = c++
FC        = f90
CFLAGS    = $(DEFINES) -pvctl,fullmsg,loopcnt=20000 -Onooverlap \
           -hfloat0 -hnostkchk
FFLAGS    = $(DEFINES) -Wf'-pvctl,fullmsg -ptr byte -b' -float0 -ew
CXXFLAGS  = -V $(DEFINES) -pvctl,loopcnt=20000,fullmsg -acct $(OPT)
RANLIB    =
CPPFLAGS += -I/usr/include/CC
OPT       = -O nooverlap -Nstkchk
NOVEC     = -c -Nstkchk -acct -dir novec
FORTLIBS  = -llapack_64 -lblas_64 -lfft_64
ifdef MPI
  LIBS    = -lfem_mp -lalp $(FORTLIBS) -lm -lC -lmpi -lpthread
else
  LIBS    = -lfem -lalp $(FORTLIBS) -lm -lC
endif
LDFLAGS   =  -verbose -V -p -f90libew $(CXXFLAGS) $(CPPFLAGS) $(LIBDIR) $(LIBS)
CLDFLAGS  = $(LDFLAGS)
endif

# ----------------------------------------------------------------------------
# Installation of header files.
#
default:
	rm -f   ../include/Sem.h		\
		../include/Feml.h		\
		../include/Mesh.h		\
		../include/Geometry.h 
	cp Sem.h Feml.h Mesh.h Geometry.h ../include

#-----------------------------------------------------------------------------
# Force installation of headers from another make.
#
sem:
	cd $(VPATH); gmake default

# -----------------------------------------------------------------------------
# Clean up.
#
clean:
	rm -f *.o *.O *.a *~ *.ti SunWS_cache
	rm -rf ILDUMPS ii_files ptrepository
