##############################################################################
# Makefile for lns, spectral element solver for linearised unsteady
# Navier--Stokes problems, and arnoldi, a Krylov-method solver for
# stability problems (including Floquet problems).
#
# $Id$
##############################################################################

SEM   = ../semtex
VPATH = $(SEM)/src

ifdef FLIP		# -- Use the "half-period-flip map".
  KPROG = arnoldi-H
  DEFINES += -DFLIP
else
  KPROG = arnoldi
endif

ifdef TBCS		# -- Recompute mode 0 velocity BCs every time step.
  DEFINES += -DTBCS
  KPROG := $(KPROG)_tbcs
endif

ifdef OLDCODE		# -- Make the old cylindrical code version.
  DEFINES += -DOLDCODE
  KPROG   := $(KPROG)_old
endif

# ----------------------------------------------------------------------------
# Semtex standard sources.
#
SEMFILES = analysis auxfield BCmgr boundary bsys condition domain	\
	   edge element family feml field geometry history integration	\
	   matrix mesh misc particle pressure
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = sem_h

# ----------------------------------------------------------------------------
# Build eigensolver code.

FL_FILES = arnoldi integrate stabanalysis
FL_OBJ   = $(addsuffix .o,$(FL_FILES))

$(KPROG): sem $(FL_OBJ) $(SEMOBJ) stab_h
	$(LD) -o $@ $(FL_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build linearised Navier-Stokes solver.

ST_FILES = lns integrate stabanalysis
ST_OBJ   = $(addsuffix .o,$(ST_FILES))

lns: sem $(ST_OBJ) $(SEMOBJ) stab_h
	$(LD) -o $@ $(ST_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build combine, a utility to combine a base and a perturbation field.
#
combine: combine.o
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build normalise, a utility to normalise eigenvectors (2-norm based).
#
normalise: normalise.o
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)


# ----------------------------------------------------------------------------
# Build flipmap, a utility to find a symmetry-related mapping of the mesh.
#
flipmap: flipmap.C
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

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
