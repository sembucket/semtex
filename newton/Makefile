##############################################################################
# Makefile for newton, spectral element solver for steady
# Navier--Stokes problems.
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

N_FILES = newton advect integrate bicgstab
N_OBJ   = $(addsuffix .o,$(N_FILES))

newton: sem $(N_OBJ) $(SEMOBJ) newt.h
	$(LD) -o $@ $(N_OBJ) $(SEMOBJ) $(LDFLAGS)

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
