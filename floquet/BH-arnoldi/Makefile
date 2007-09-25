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

ifdef ARPACK
  DEFINES += -DARPACK
endif

CPPFLAGS = -I. -I$(SEM)/include

# ----------------------------------------------------------------------------
# Build arnoldi solver:

OBJ = iohb.o dspmvc.o arnoldi.o


arnoldi:	$(OBJ)
	$(LD) -o $@ $(OBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

ifdef ARPACK
  F77LIBS = -larpack -framework Accelerate -lg95
endif

