##############################################################################
# Makefile for aero, spectral element solver for unsteady Navier--Stokes
# problems in accelerating reference frame.
#
# $Id$
##############################################################################

ifdef MPI
  PROG = aero_mp
else
  PROG = aero
endif

ifdef LES		# -- Compile in Smagorinsky-based LES.
  DEFINES += -DLES
endif

ifdef DAMPING		# -- Compile in van-Driest damping (caveat emptor).
  DEFINES += -DDAMPING
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
# Default build for aeroelastic Navier-Stokes solver.

NS_OBJ = aeroanalysis.o body.o drive.o NS.o 
NS_HDR = aero.h
ifdef LES
  NS_OBJ += eddyvis.o
endif


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

# ----------------------------------------------------------------------------
# Refreq is a utility to compute motion phase angles at a given time.

refreq: refreq.C
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Dyres is a utility to add body-motion variables to the end of a .bdy file.

dyres: dyres.c
	$(CC) -o $@ $@.c

# ----------------------------------------------------------------------------
# Relv is a utility to to add frame velocity to velocities in field dumps.

relv: relv.c
	$(CC) -o $@ $@.c $(CPPFLAGS) $(LDFLAGS)
