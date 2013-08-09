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

ifdef STOKES
  DEFINES += -DSTOKES
endif

ifdef ALIAS		# -- No dealiasing done.
  DEFINES += -DALIAS
endif

# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield bcmgr boundary bsys condition domain	\
	    edge element family feml field geometry history integration \
	    matrix mesh misc particle pressure statistics
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = sem.h

# ----------------------------------------------------------------------------
# Default build.

NS_OBJ = scatanalysis.o drive.o integrate.o nonlinear.o integrates.o
NS_HDR = scat.h

$(PROG): sem $(SEMOBJ) $(NS_OBJ) 
	$(LD) -o $@ $(NS_OBJ) $(SEMOBJ) $(LDFLAGS)

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Specialised version of rstress

rstress: rstress.c
	$(CC) -o $@_s rstress.c $(CFLAGS) $(CPPFLAGS) $(CLDFLAGS) -lm

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
