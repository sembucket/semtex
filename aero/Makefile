##############################################################################
# Makefile for aero, spectral element solver for unsteady Navier--Stokes
# problems in accelerating reference frame.
#
# $Id$
##############################################################################

SEM      = ..
VPATH    = $(SEM)/src

ARCH     = $(shell uname -s)

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
endif

CPPFLAGS = -I. -I$(SEM)/include

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x
CXX	 = CC
LD	 = CC
CC	 = cc
CXXFLAGS = $(DEFINES) -fullwarn -O2
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -lfpe -L$(SEM)/lib/$(ARCH) -lfem -lalp \
	   $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CXX	 = CC
LD	 = CC
CC       = cc
CXXFLAGS =  $(DEFINES) -O2
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS)  -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 AlphaServer multiprocessor
CXX	 = cxx
LD	 = cxx
CC       = cc
CXXFLAGS =  $(DEFINES) -migrate -tune host -O2
F77LIBS  = -ldxml
LDFLAGS  = $(CXXFLAGS) -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CXX	 = g++
LD	 = g++
CC       = gcc
CXXFLAGS = $(DEFINES) -O3
F77LIBS  = -llapack -lblas -lf2c
LDFLAGS  = -gstabs+ -L$(SEM)/lib/$(ARCH) -L/usr/local/lib	\
	   -lfem -lalp $(F77LIBS)
endif

ifeq ($(ARCH),UNIX_System_V)    # -- Fujitsu VPP300.
LIBDIR   = -L$(SEM)/lib/$(ARCH) 
LIBDIR  += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
LIBS     = -lfem -lalp
LIBS    += -llapack -lblasvpp  -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
CC       = vcc
CXX      = /usr/uxp/C++/bin/CC
LD       = /usr/uxp/C++/bin/CC
#OPTIM    = -Wv,-m3,-Ps 
OPTIM    = -Wv,-m3,-Om,-te
CXXFLAGS = $(DEFINES)  $(OPTIM)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS)
endif

ifeq ($(ARCH),SUPER-UX)         # -- NEC SX-4.
LIBDIR   = -L$(SEM)/lib/$(ARCH) 
LIBDIR  += 
LIBS     = -lfem -lalp -llapack -lblas	\
           -lu77sxe -lp77sx -lv77sxe -lf77sxe -li77sx -lm
CC       = cc
FC       = f77
CXX      = CC
LD       = CC
OPTIM    = -g -hacct
CXXFLAGS = $(DEFINES)  $(OPTIM)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS) -hacct
endif

# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield BCmgr boundary condition domain	\
	    element feml field geometry history integration 	\
	    matrix mesh misc particle pressure statistics
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h

# ----------------------------------------------------------------------------
# Default build for aeroelastic Navier-Stokes solver.

NS_OBJ = aeroanalysis.o body.o drive.o NS.o 
NS_HDR = aero.h

aero:	sem $(SEMOBJ) $(NS_OBJ) 
	$(LD) -o $@ $(NS_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
#
$(SEMOBJ): $(SEMHDR)

#-----------------------------------------------------------------------------
sem:
	cd $(VPATH); gmake install

# ----------------------------------------------------------------------------
# Refreq is a utility to compute motion phase angles at a given time.

refreq: refreq.C
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Dyres is a utility to add frame velocity to velocities in field dumps.

dyres: dyres.c
	$(CC) -o $@ $@.c

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
	rm -rf ptrepository
