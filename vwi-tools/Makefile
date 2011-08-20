##############################################################################
# Makefile for VWI utilities.
#
# $Id$
##############################################################################

SEM   = ../semtex
VPATH = $(SEM)/src

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile
SEMHDR = sem.h

# ----------------------------------------------------------------------------
# forc makes 2D roll forcing from streak and wave solutions.
#

FOBJ = feml.o mesh.o element.o family.o auxfield.o geometry.o data2df.o

forc: forc.o $(FOBJ)
	$(CXX) $(CPPFLAGS) -o $@ $@.o $(FOBJ) $(LDFLAGS)

# -----------------------------------------------------------------------------
# wavestress makes 2C Reynolds stress components from a 2C complex mode.
#

wavestress: wavestress.c
	$(CC) $(CPPFLAGS) -o $@ $@.c $(LDFLAGS)

# -----------------------------------------------------------------------------
# deltaweight multiplies delta function weight approximation onto roll forcing.
#

deltaweight: deltaweight.c
	$(CC) $(CPPFLAGS) -o $@ $@.c $(LDFLAGS)
