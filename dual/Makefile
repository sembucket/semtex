##############################################################################
# Makefile for dual, spectral element solver for unsteady
# Navier--Stokes problems with 2 Fourier modes.  Single-processor
# only.
#
# $Id$
##############################################################################

SEM   = ../semtex
VPATH = $(SEM)/src
PROG  = dual

# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield BCmgr boundary bsys condition domain \
	    element feml field geometry history integration  	   \
	    matrix mesh misc particle pressure statistics
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
#
SOLOBJ = drive.o NS.o

$(PROG): sem $(SEMOBJ) $(SOLOBJ) 
	$(LD) -o $@ $(SOLOBJ) $(SEMOBJ) $(LDFLAGS)

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Dual2sem converts dual field files to sem format & vice-versa.
#
dual2sem: dual2sem.C
	$(CXX) $(CPPFLAGS) -g -o $@ $@.C $(LDFLAGS)

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
