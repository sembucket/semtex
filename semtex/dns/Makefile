##############################################################################
# Makefile for dns, spectral element solver for unsteady Navier--Stokes
# problems.
#
# $Id$
##############################################################################

PROG = dns

SEM   = ..
VPATH = $(SEM)/src

ifdef STOKES		# -- Turn off nonlinear terms; create Stokes solver.
  DEFINES += -DSTOKES
  PROG := $(PROG)_stokes
endif

ifdef ALIAS		# -- No dealiasing (parallel is never dealiased).
  DEFINES += -DALIAS
  PROG := $(PROG)_alias
endif

ifdef SKEW		# -- Make skew-symmetric nonlinear terms.
  DEFINES += -DSKEW
  PROG := $(PROG)_skew
endif

ifdef TBCS		# -- Recompute mode 0 velocity BCs every time step.
  DEFINES += -DTBCS
  PROG := $(PROG)_tbcs
endif

ifdef MPI		# -- Make a parallel executable.
  PROG := $(PROG)_mp
endif

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
NS_FILES = drive integrate dnsanalysis
NS_OBJ   = $(addsuffix .o,$(NS_FILES))

$(PROG): sem $(SEMOBJ) $(NS_OBJ) 
	$(LD) -o $@ $(NS_OBJ) $(SEMOBJ) $(LDFLAGS)

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
