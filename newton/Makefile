##############################################################################
# Makefile for newton, spectral element solver for steady
# Navier--Stokes problems.
#
# $Id$
##############################################################################

SEM   = ../semtex
VPATH = $(SEM)/src

ifdef NSPCG
  F77LIBS += -lnspcg
  DEFINES += -DNSPCG
endif

# ----------------------------------------------------------------------------
# Semtex standard sources.
#
SEMFILES = auxfield BCmgr boundary bsys condition domain	\
	   element feml field geometry history integration	\
	   matrix mesh misc particle pressure statistics
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = Sem.h

# ----------------------------------------------------------------------------
# Build linearised Navier-Stokes solver.

N_FILES = newton advect integrate bicgstab
N_OBJ   = $(addsuffix .o,$(N_FILES))
#	$(LD) -o $@ $(N_OBJ) $(SEMOBJ) $(LDFLAGS) -L. -lnspcg -lftn

newton: sem $(N_OBJ) $(SEMOBJ) newt.h
	$(LD) -o $@ $(N_OBJ) $(SEMOBJ) -L. -lnspcg -L../semtex/lib/IRIX64 -lfem -lalp -lcomplib.sgimath -lftn -lm 

# ----------------------------------------------------------------------------
# Header dependency.

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
