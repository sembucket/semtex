##############################################################################
# Makefile for scat, spectral element solver for unsteady
# Navier--Stokes problems, with passive SCAlar Transport.
#
# $Id$
##############################################################################

ifdef MPI
  PROG = scat_mp
else
  PROG = scat
endif

SEM   = ../semtex
VPATH = $(SEM)/src

# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield BCmgr boundary bsys condition domain	\
	    element feml field geometry history integration 		\
	    matrix mesh misc particle pressure statistics
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h

# ----------------------------------------------------------------------------
# Default build.

NS_OBJ = scatanalysis.o drive.o integrate.o
NS_HDR = scat.h

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
