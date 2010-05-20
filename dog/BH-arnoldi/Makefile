##############################################################################
# Makefile for arnoldi eigensolver.
#
# $Id$
##############################################################################

SEM      = ../../semtex
VPATH    = $(SEM)/src

ARCH     = $(shell uname -s)

DEFINES = 

ifdef DEBUG
  DEFINES += -DDEBUG
endif

PROG = arnoldi-DB

ifdef ARPACK
  DEFINES += -DARPACK
  PROG = arnoldi-AR
endif

CPPFLAGS = -I. -I$(SEM)/include

# ----------------------------------------------------------------------------
# Build arnoldi solver:

OBJ = iohb.o dspmvc.o arnoldi.o

$(PROG):	$(OBJ)
	$(LD) -o $@ $(OBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

ifdef ARPACK
#  F77LIBS = -larpack -llapack -lfastblas -lgfortran
  F77LIBS := -larpack $(F77LIBS)
endif

