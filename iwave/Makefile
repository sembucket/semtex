##############################################################################
# Makefile for dns, spectral element solver for unsteady Navier--Stokes
# problems.
#
# $Id$
##############################################################################

ifdef MPI
  PROG = iwave_mp
else
  PROG = iwave
endif

SEM   = ../semtex
VPATH = $(SEM)/src

ifdef STOKES
  DEFINES += -DSTOKES
endif

# ----------------------------------------------------------------------------
#
SEMFILES = analysis auxfield BCmgr boundary bsys condition domain   \
	   element feml field geometry history integration  	\
	   matrix mesh misc particle pressure statistics
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = Sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
NS_OBJ = drive.o NS.o dnsanalysis.o

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
