##############################################################################
# Makefile for dual, spectral element solver for unsteady
# Navier--Stokes problems with 2 Fourier modes.  Single-processor
# only.
#
# $Id$
##############################################################################

PROG  = dual

SEM   = ../semtex
VPATH = $(SEM)/src

ifdef STOKES
  DEFINES += -DSTOKES
  PROG    := $(PROG)_stokes
endif

ifdef SKEW
  DEFINES += -DSKEW
  PROG    := $(PROG)_skew
endif

#
# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield bcmgr boundary bsys condition domain	\
            edge element family feml field geometry integration		\
	    matrix mesh misc pressure statistics
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
#
SOLOBJ = drive.o NS.o dualanalysis.o

$(PROG): sem $(SOLOBJ) $(SEMOBJ)
	$(LD) -o $@ $(SOLOBJ) $(SEMOBJ) $(LDFLAGS)

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Dual2sem converts dual field files to sem format & vice-versa.
#
dual2sem: dual2sem.o
	$(CXX) $(CPPFLAGS) -o $@ $@.o $(LDFLAGS)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

#-----------------------------------------------------------------------------
# Machine-dependent compilations.
#
ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif
