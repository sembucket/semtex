##############################################################################
# Makefile for dns, spectral element solver for unsteady Navier--Stokes
# problems.
#
# $Id$
##############################################################################

SEM   = ../semtex
VPATH = $(SEM)/src

DEFINES += -DSTABILITY

# ----------------------------------------------------------------------------
# Semtex standard sources.
#
SEMFILES = analysis auxfield BCmgr boundary bsys condition domain	\
	   element feml field geometry history integration		\
	   matrix mesh misc particle pressure statistics
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = Sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
#

ST_FILES = stab NS stabanalysis
ST_OBJ   = $(addsuffix .o,$(ST_FILES))

FL_FILES = arnoldi NS stabanalysis
FL_OBJ   = $(addsuffix .o,$(FL_FILES))

stab: sem $(ST_OBJ) $(SEMOBJ)
	$(LD) -o $@ $(ST_OBJ) $(SEMOBJ) $(LDFLAGS)

flok: sem $(FL_OBJ) $(SEMOBJ)
	$(LD) -o $@ $(FL_OBJ) $(SEMOBJ) $(LDFLAGS)

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
