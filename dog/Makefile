##############################################################################
# Makefile for dns, spectral element solver for unsteady Navier--Stokes
# problems.
#
# $Id$
##############################################################################

SEM   = ../semtex
VPATH = $(SEM)/src

# ----------------------------------------------------------------------------
# Semtex standard sources.
#
SEMFILES = analysis auxfield BCmgr boundary bsys condition domain	\
	   element feml field geometry history integration		\
	   matrix mesh misc particle pressure statistics
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = Sem.h

# ----------------------------------------------------------------------------
# Build linearised Navier-Stokes solver.

ST_FILES = lns integrate stabanalysis
ST_OBJ   = $(addsuffix .o,$(ST_FILES))

lns: sem $(ST_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(ST_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build eigensolver code

FL_FILES = arnoldi integrate stabanalysis
FL_OBJ   = $(addsuffix .o,$(FL_FILES))

arnoldi: sem $(FL_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(FL_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Header dependency.

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
