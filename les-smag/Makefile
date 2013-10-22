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

ifdef NOMODEL		# -- This will produce dns-like execution.
  DEFINES += -DNOMODEL
endif

ifdef DAMPING		# -- Compile in van-Driest damping (caveat emptor).
  DEFINES += -DDAMPING
endif

ifdef ALIAS		# -- No dealiasing done.
  DEFINES += -DALIAS
endif

# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield bcmgr boundary bsys condition domain	\
	    edge element family feml field geometry history integration	\
	    matrix mesh misc particle pressure statistics svv
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = sem.h

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
