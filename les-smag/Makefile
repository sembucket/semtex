##############################################################################
# Makefile for les, spectral element solver for unsteady
# Navier--Stokes problems incorporating LES in Cartesian geometries.
#
# $Id$
##############################################################################

SEM   = ../semtex
VPATH = $(SEM)/src

ifdef MPI
  PROG = les-smag_mp
else
  PROG = les-smag
endif

ifdef DAMPING		# -- Compile in van-Driest damping (caveat emptor).
  DEFINES += -DDAMPING
endif

ifdef ALIAS		# -- No dealiasing done.
  DEFINES += -DALIAS
endif

# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield BCmgr boundary bsys condition domain \
	    element feml field geometry history integration  	   \
	    matrix mesh misc particle pressure statistics
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
LESOBJ = drive.o NS.o eddyvis.o lesanalysis.o

$(PROG): sem $(SEMOBJ) $(LESOBJ) 
	$(LD) -o $@ $(LESOBJ) $(SEMOBJ) $(LDFLAGS)

$(SEMOBJ): $(SEMHDR)

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
