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

# ----------------------------------------------------------------------------
# Change variables to use NSPCG library, assumed linked to local directory.
#

ifdef NSPCG
  F77LIBS := -L. -lnspcg $(F77LIBS)
  DEFINES += -DNSPCG
endif

#-----------------------------------------------------------------------------
# Machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif
