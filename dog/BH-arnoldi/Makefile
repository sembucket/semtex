##############################################################################
# Makefile for arnoldi eigensolver.
#
# $Id$
##############################################################################

SEM      = ../../semtex
VPATH    = $(SEM)/src

ARCH     = $(shell uname -s)

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
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
