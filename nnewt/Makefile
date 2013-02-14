##############################################################################
# Makefile for nnewt, spectral element solver for DNS of unsteady
# Navier--Stokes problems incorporating non-Newtonian fluids.
#
# $Id$
##############################################################################

PROG = nnewt

SEM   = ../semtex
VPATH = $(SEM)/src

ifdef MPI
  PROG := $(PROG)_mp
endif

ifdef NOMODEL
  DEFINES += -DNOMODEL
  PROG    := $(PROG)_nomodel
endif

# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield bcmgr boundary bsys condition domain   \
	    edge element family feml field geometry history \
	    integration matrix mesh misc particle pressure statistics
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
NNEWTOBJ = drive.o NS.o viscosity.o nnewtanalysis.o

$(PROG): sem $(SEMOBJ) $(NNEWTOBJ)  nnewt.h
	$(LD) -o $@ $(NNEWTOBJ) $(SEMOBJ) $(LDFLAGS)

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Build eneq_nn, which computes terms in fluctuating flow energy equations.
#
ENOBJ = eneq.o feml.o mesh.o element.o family.o auxfield.o geometry.o misc.o
$(ENOBJ): $(SEMHDR)

eneq: $(ENOBJ)
	$(CXX) $(CPPFLAGS) -o $@_nn $(ENOBJ) $(LDFLAGS)

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
