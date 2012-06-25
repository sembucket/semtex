##############################################################################
# Makefile for optimal BC calculations, based on Xuerui Mao's
# modifiactions to dog, 2012.
#
# $Id$
##############################################################################
SEM   = ../semtex
UTIL  = $(SEM)/utility

# where to find the files not present in the local directory.
VPATH = $(SEM)/src

ifdef TBCS		# -- Recompute mode 0 velocity BCs every time step.
  DEFINES += -DTBCS
  KPROG := $(KPROG)_tbcs
endif

ifdef ARPACK
  DEFINES += -DARPACK
  KPROG   := $(KPROG)-AR
endif

# ----------------------------------------------------------------------------
# Semtex standard sources.
# 
SEMFILES = analysis auxfield bcmgr boundary bsys condition data2df	\
	   domain edge element family feml field geometry history	\
           integration matrix mesh misc particle pressure
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = $(addsuffix .h,$(SEMFILES)) sem.h

# ----------------------------------------------------------------------------
# Build opt solver for boundary perturbations.

FL_OPTBC_FILES = optbc integrate stabanalysis
FL_OPTBC_OBJ   = $(addsuffix .o,$(FL_OPTBC_FILES))

optbc: sem $(FL_OPTBC_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(FL_OPTBC_OBJ) $(SEMOBJ) $(LDFLAGS)
	rm $(UTIL)/enumerate.o $(UTIL)/enumerate;
	cp mesh.o $(UTIL)/mesh.o; cp feml.o $(UTIL)/feml.o;
	cd $(UTIL); $(MAKE) enumerate

# ----------------------------------------------------------------------------
# minimize the unsteadiness such as vortex shedding.

FL_UNSTEADY_FILES = unsteady integrate_unsteady stabanalysis
FL_UNSTEADY_OBJ   = $(addsuffix .o,$(FL_UNSTEADY_FILES))

unsteady: sem $(FL_UNSTEADY_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(FL_UNSTEADY_OBJ) $(SEMOBJ) $(LDFLAGS)
# ----------------------------------------------------------------------------
# Build bc solver to minimize energy of ic.

FL_ICBC_FILES = icbc integrate stabanalysis
FL_ICBC_OBJ   = $(addsuffix .o,$(FL_ICBC_FILES))

icbc: sem $(FL_ICBC_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(FL_ICBC_OBJ) $(SEMOBJ) $(LDFLAGS)
# ----------------------------------------------------------------------------
# Build opt solver for initial perturbations.

FL_OPTIC_FILES = optic integrate stabanalysis
FL_OPTIC_OBJ   = $(addsuffix .o,$(FL_OPTIC_FILES))

optic: sem $(FL_OPTIC_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(FL_OPTIC_OBJ) $(SEMOBJ) $(LDFLAGS)
# ----------------------------------------------------------------------------
# Build eig solver for initial perturbations.

FL_EIGIC_FILES = eigic integrate stabanalysis
FL_EIGIC_OBJ   = $(addsuffix .o,$(FL_EIGIC_FILES))

eigic: sem $(FL_EIGIC_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(FL_EIGIC_OBJ) $(SEMOBJ) $(LDFLAGS)
# ----------------------------------------------------------------------------
# Build controllability code

FL_CABILITY_FILES = cability integrate stabanalysis
FL_CABILITY_OBJ   = $(addsuffix .o,$(FL_CABILITY_FILES))

cability: sem $(FL_CABILITY_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(FL_CABILITY_OBJ) $(SEMOBJ) $(LDFLAGS)
# ----------------------------------------------------------------------------
# Build eigenvalue solver for boundary perturbations.

FL_EIGBC_FILES = eigbc integrate stabanalysis
FL_EIGBC_OBJ   = $(addsuffix .o,$(FL_EIGBC_FILES))

eigbc: sem $(FL_EIGBC_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(FL_EIGBC_OBJ) $(SEMOBJ) $(LDFLAGS)
# ----------------------------------------------------------------------------
# Build linearised Navier-Stokes solver.

ST_FILES = lnsbc integrate stabanalysis
ST_OBJ   = $(addsuffix .o,$(ST_FILES))

lnsbc: sem $(ST_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(ST_OBJ) $(SEMOBJ) $(LDFLAGS)

#----------------------------------------------------------------------------
# Build diffbcs, calculate the difference between two boundary perturbations.

ST_DIFFBCS_FILES = diffbcs  stabanalysis
ST_DIFFBCS_OBJ   = $(addsuffix .o,$(ST_DIFFBCS_FILES))

diffbcs: sem $(ST_DIFFBCS_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(ST_DIFFBCS_OBJ) $(SEMOBJ) $(LDFLAGS)

#----------------------------------------------------------------------------
# Build diffics, calculate the difference between two initial perturbations.

ST_DIFFICS_FILES = diffics  stabanalysis
ST_DIFFICS_OBJ   = $(addsuffix .o,$(ST_DIFFICS_FILES))

diffics: sem $(ST_DIFFICS_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(ST_DIFFICS_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Header dependency.

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

ifdef ARPACK
  F77LIBS := -larpack $(F77LIBS)
endif

#-----------------------------------------------------------------------------
# Machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif
