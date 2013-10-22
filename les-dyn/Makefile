##############################################################################
# Makefile for les-dyn, spectral element solver for unsteady
# Navier--Stokes problems incorporating dynamic LES in Cartesian
# geometries.
#
# $Id$
##############################################################################

SEM   = ../semtex
VPATH = $(SEM)/src

# ----------------------------------------------------------------------------
# Conditional compilation flags.

ifdef MPI			# -- Link parallel code.
  PROG = les-dyn_mp
else
  PROG = les-dyn
endif

ifdef NOMODEL			# -- This will produce DNS.
  DEFINES += -DNOMODEL
endif

ifdef SMAG			# -- Over-ride dynamic estimate with C_SMAG.
  DEFINES += -DSMAG
endif

ifdef PROJ			# -- Use 2D projection instead of BVD filter.
  DEFINES += -DPROJ
endif

# ----------------------------------------------------------------------------
# Standard semtex objects.

SEMFILES = analysis auxfield bcmgr boundary bsys condition domain \
	   edge element family feml field geometry history integration \
	   matrix mesh misc particle pressure statistics svv
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = sem.h

# ----------------------------------------------------------------------------
# Local objects and build.

LESFILES = drive integrate nonlinear filter lesanalysis sumintegrator
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
