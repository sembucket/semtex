##############################################################################
# Makefile for chroma, specialised spectral element solver for unsteady
# 3D advection--diffusion of scalar in WCOT chromatographs.
#
# $Id$
##############################################################################

SEM   = ../semtex
VPATH = $(SEM)/src

ARCH = $(shell uname -s)

ifeq ($(ARCH),Linux)
  MACH = $(shell uname -m)
  ifeq ($(MACH),alpha)
    ARCH = Linux-alpha
  endif
  ifeq ($(MACH),i586)
    ARCH = Linux-i586
  endif
endif

ifdef DEBUG
  DEFINES = -DDEBUG
else
  DEFINES  =
endif

ifdef MPI
  PROG    = chroma_mp
else
  PROG = chroma
endif

# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield BCmgr boundary bsys condition domain   \
	    element feml field geometry history integration \
	    matrix mesh misc pressure
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h

# ----------------------------------------------------------------------------
# Default build for solver.
NS_OBJ = drive.o patch.o transport.o

$(PROG): sem $(SEMOBJ) $(NS_OBJ) 
	$(LD) -o $@ $(NS_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build post-processor for analysis of integrated concentration profiles.

momenta: sem momenta.C
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)


# ----------------------------------------------------------------------------
#
$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

#-----------------------------------------------------------------------------
# Machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif
