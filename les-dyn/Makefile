##############################################################################
# Makefile for les-dyn, spectral element solver for unsteady
# Navier--Stokes problems incorporating dynamic LES in Cartesian
# geometries.
#
# $Id$
##############################################################################

ifdef MPI
  PROG = les-dyn_mp
else
  PROG = les-dyn
endif

SEM   = ../semtex
VPATH = $(SEM)/src

ifdef NOMODEL			# -- This will produce DNS.
  DEFINES += -DNOMODEL
endif

ifdef SMAG			# -- Over-ride dynamic estimate with C_SMAG
  DEFINES += -DSMAG
endif

# ----------------------------------------------------------------------------
#
SEMFILES = analysis auxfield BCmgr boundary bsys condition domain \
	   element feml field geometry history integration  	  \
	   matrix mesh misc particle pressure statistics
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = Sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
LESFILES = drive integrate nonlinear filter lesanalysis
LESOBJ   = $(addsuffix .o,$(LESFILES))

$(PROG): sem $(SEMOBJ) $(LESOBJ) 
	$(LD) -o $@ $(LESOBJ) $(SEMOBJ) $(LDFLAGS)

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
