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

ifdef STOKES
  DEFINES += -DSTOKES
endif

#
# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield BCmgr boundary bsys condition	\
	    domain element feml field geometry integration	\
	    matrix mesh misc pressure statistics
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
#
SOLOBJ = drive.o NS.o dualanalysis.o

$(PROG): sem $(SEMOBJ) $(SOLOBJ) 
	$(LD) -o $@ $(SOLOBJ) $(SEMOBJ) $(LDFLAGS)

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Dual2sem converts dual field files to sem format & vice-versa.
#
dual2sem: dual2sem.o
	$(CXX) $(CPPFLAGS) -o $@ $@.o $(LDFLAGS)

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
