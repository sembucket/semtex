##############################################################################
# Makefile for aero, spectral element solver for unsteady Navier--Stokes
# problems in accelerating reference frame.
#
# $Id$
##############################################################################

FEM      = ..
VPATH    = $(FEM)/src

ARCH     = $(shell uname -s)

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
endif

CPPFLAGS = -I. -I$(FEM)/include

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -fullwarn -g
F77LIBS	 = -lcomplib.sgimath -lisam
LDFLAGS  = -lfpe -L$(FEM)/lib/$(ARCH) -lfem -lalp \
	   $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CXX	 = CC
LD	 = CC
CXXFLAGS =  $(DEFINES) -O3
F77LIBS	 = -lcomplib.sgimath -lisam
LDFLAGS  = -L$(FEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 AlphaServer multiprocessor
CXX	 = cxx
LD	 = cxx
CXXFLAGS =  $(DEFINES) -migrate -tune host -O2
F77LIBS  = -ldxml
LDFLAGS  = $(CXXFLAGS) -L$(FEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CXX	 = g++
LD	 = g++
CXXFLAGS = $(DEFINES) -O3
F77LIBS  = -llapack -lblas -lf2c
LDFLAGS  = -gstabs+ -L$(SEM)/lib/$(ARCH) -L/usr/local/lib	\
	   -lfem -lalp $(F77LIBS)
endif

# ----------------------------------------------------------------------------
#
SEMFILES  = auxfield BCmgr boundary condition domain element	\
	    feml field integration matrix mesh misc pressure
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h

# ----------------------------------------------------------------------------
# Default build for aeroelastic Navier-Stokes solver.

NS_OBJ = drive.o NS.o analysis.o body.o
NS_HDR = aero.h

aero:	$(SEMOBJ) $(NS_OBJ) 
	cd ../src; $(MAKE)
	$(LD) -o $@ $(NS_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
#
$(SEMOBJ): $(SEMHDR)
$(NS_OBJ): $(NS_HDR)

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
