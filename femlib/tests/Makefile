##############################################################################
# Compile test codes for femlib library functions.
#
# $Id$
##############################################################################

SEMHOME  = ../..
LIB      = libfem.a

ARCH = $(shell uname -s)
MACH = $(shell uname -m)

ifeq ($(ARCH),Linux)
  ifeq ($(MACH),alpha)
    ARCH = Linux-alpha
  endif
  ifeq ($(MACH),i586)
    ARCH = Linux-i586
  endif
  ifeq ($(MACH),i686)
    ARCH = Linux-i686
  endif
endif

DEFINES  =
ifdef DEBUG
  DEFINES += -DDEBUG
endif

CC	 = cc
CXX	 = CC
FC	 = f77
OPTIONS  = -g
CFLAGS	 = $(OPTIONS)
FFLAGS   = $(OPTIONS)
ARFLAGS	 = cr
CPPFLAGS = -I$(SEMHOME)/include
RANLIB   =

# -------------- Compiler options for specific machines. ---------------------

ifeq ($(ARCH),IRIX)			# -- SGI IRIX 32-bit
CFLAGS     = $(DEFINES) -fullwarn -O2
FFLAGS     = $(DEFINES) -O2
endif

ifeq ($(ARCH),IRIX64)			# -- SGI IRIX 64-bit
CXXFLAGS   = $(DEFINES) -mips4 -64 -g 
CFLAGS     = $(DEFINES) -mips4 -64 -g 
FFLAGS     = $(DEFINES) -mips4 -64 -g
LDFLAGS    = $(CXXFLAGS)
F77LIBS    = -lcomplib.sgimath -lftn
endif

ifeq ($(ARCH),OSF1)			# -- DEC Alpha OSF1
FC	   = f90
CXX	   = cxx
CC         = cc
CFLAGS     = $(DEFINES) -fast -inline speed -tune host -O4
FFLAGS     = $(DEFINES) -fast -inline speed -tune host -O5 -unroll 2
ifdef MPI
  LIB       = libfem_mp.a
  DEFINES  += -DMPI
  CPPFLAGS += -I/usr/local/mpi/include
  CFLAGS   += -pthread
  FFLAGS   += -D_THREAD_SAFE
endif
endif

ifeq ($(ARCH),Linux-i686)		# -- Linux Intel system
CC         = gcc
FC         = g77
CXX	   = g++
CFLAGS     = $(DEFINES) -O3
FFLAGS     = $(DEFINES) -O3
F77LIBS  = -llapack -lblas -lg2c
endif

ifeq ($(ARCH),Linux-alpha)		# -- Linux alpha system.
CC	 = ccc
CXX	 = g++
LD	 = g++
ifdef DEBUG
  OPT    = -g
else
  OPT    = -O2
endif
CXXFLAGS = $(DEFINES) $(OPT)
F77LIBS  = -lcxml -lots -lfor -lg2c
LDFLAGS  = -L$(SEM)/lib/$(ARCH) -L/usr/local/lib	\
	   -lfem -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),UNIX_System_V)            # -- Fujitsu VPP300
CC         = vcc
FC         = frt
CFLAGS     = $(DEFINES) -Wv,-m3,-Om,-te
FFLAGS     = $(DEFINES) -Wv,-m3,-te
ifdef MPI
  LIB       = libfem_mp.a
  DEFINES  += -DMPI
  CPPFLAGS += -I/usr/lang/mpi/include
endif
endif

ifeq ($(ARCH),SUPER-UX)                 # -- NEC SX-4
CFLAGS     = $(DEFINES) -pvctl,fullmsg -hfloat0 -hnostkchk
FFLAGS     = $(DEFINES) -Cvopt -Wf'-pvctl,fullmsg -ptr byte -b' -float0 -ew
ifdef MPI
  LIB       = libfem_mp.a
  DEFINES  += -DMPI
  CFLAGS   += 
  FFLAGS   += -P multi
endif
endif

testpoly: testpoly.c
	cc $(CFLAGS) -o testpoly testpoly.c -I$(SEMHOME)/include -I. \
	-L$(SEMHOME)/lib -lfem -lalp -lm

testFFT: testFFT.C
	$(CXX) -o $@ testFFT.C -I$(SEMHOME)/include -I. -g \
	$(CXXFLAGS) -L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lftn -lm

testDLT: testDLT.C
	$(CXX) -o $@ testDLT.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp -lm

testDLT2: testDLT2.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

testDLT3: testDLT3.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
