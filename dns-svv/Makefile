##############################################################################
# Makefile for dns, spectral element solver for unsteady Navier--Stokes
# problems.
#
# $Id$
##############################################################################

PROG = dns_svv

SEM   = ../semtex
VPATH = $(SEM)/src

ifdef STOKES		# -- Turn off nonlinear terms; create Stokes solver.
  DEFINES += -DSTOKES
  PROG    := $(PROG)_stokes
endif

ifdef ALIAS		# -- No dealiasing (parallel is never dealiased).
  DEFINES += -DALIAS
  PROG    := $(PROG)_alias
endif

ifdef CONV		# -- Make convective nonlinear terms (not skew-symm).
  DEFINES += -DCONV
  PROG    := $(PROG)_conv
endif

ifdef ASS		# -- Make alternating skew-symmetric nonlinear.
  DEFINES += -DASS
  PROG    := $(PROG)_ass
endif

ifdef MPI		# -- Make a parallel executable.
  PROG := $(PROG)_mp
endif

# ----------------------------------------------------------------------------
# Semtex standard sources.
#
SEMFILES = analysis auxfield bcmgr boundary bsys condition domain   \
	   edge element family feml field geometry history \
           integration matrix mesh misc particle pressure statistics
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = $(addsuffix .h,$(SEMFILES)) sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
#
NS_FILES = drive integrate nonlinear dnsanalysis svv
NS_OBJ   = $(addsuffix .o,$(NS_FILES))

$(PROG): sem $(SEMOBJ) $(NS_OBJ) 
	$(LD) -o $@ $(NS_OBJ) $(SEMOBJ) $(LDFLAGS)

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

# ----------------------------------------------------------------------------
# Clean up executables and anything else
#
distclean:
	$(MAKE) clean
	rm -f dns_svv dns_svv_tbcs dns_svv_mp dns_svv_tbcs_mp

# ----------------------------------------------------------------------------
# Make codes for time-varying BCs
#
tbcs:
	rm -f integrate.o
	$(MAKE) TBCS=1
	rm -f integrate.o

# ----------------------------------------------------------------------------
# Make unsteady Stokes solver:
#
stokes:
	rm -f nonlinear.o
	$(MAKE) STOKES=1
	rm -f nonlinear.o

#-----------------------------------------------------------------------------
# Machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif
