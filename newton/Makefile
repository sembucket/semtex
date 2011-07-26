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
SEMFILES = auxfield bcmgr boundary bsys condition domain edge		\
           element family feml field geometry history integration	\
	   matrix mesh misc particle pressure statistics
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = sem.h

ifdef NSPCG # -- Use NSPCG library for interative solver.
  N_FILES = newton advect integrate
  F77LIBS += -L. -lnspcg $(F77LIBS)
  DEFINES += -DNSPCG
else ifdef LST # -- Use Laurette Tuckerman's implementation of bcgstab.
  N_FILES = newton advect integrate bcgstab
  DEFINES += -DLST
else # -- Default: Templates' bicgstab.
  N_FILES = newton advect integrate bicgstab
endif

N_OBJ = $(addsuffix .o,$(N_FILES))

# ----------------------------------------------------------------------------
# Build linearised Navier-Stokes solver.

newton: sem $(N_OBJ) $(SEMOBJ) newt.h
	$(LD) -o $@ $(N_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Header dependency.

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#

include $(VPATH)/Makefile

ifdef NSPCG
  F77LIBS := -L. -lnspcg $(F77LIBS)
endif

#-----------------------------------------------------------------------------
# Machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif
