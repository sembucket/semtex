##############################################################################
# Makefile for grid generation codes.
#
# $Id$
##############################################################################

#ARCH     = $(shell uname -s)
#OPTIONS  = -g
#DEFINES  =

#CXXFLAGS = $(OPTIONS) $(DEFINES)
#CPPFLAGS = -I$(HOME)/include -I.

SEM=../../develop/semtex
VPATH=$(SEM)/src

# ----------------------------------------------------------------------------
# Build cgen

cgen:	cgen.o airfoil.o stretch.o rk.o odeint.o
	$(LD) -o $@ cgen.o airfoil.o stretch.o rk.o odeint.o $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build cgrid

cgrid: cgrid.f
	$(F77) -o $@ $@.f

# ----------------------------------------------------------------------------
# Build testfoil

testfoil: testfoil.o airfoil.o rk.o odeint.o
	$(LD) -o $@ $(CXXFLAGS) testfoil.o airfoil.o rk.o odeint.o $(LDFLAGS)

# ---------------------------------------------------------------------------
# Build plate

plate: plate.o stretch.o
	$(LD) -o $@ $(CXXFLAGS) plate.o stretch.o $(LDFLAGS)

# ----------------------------------------------------------------------------
# Clean up.

clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files


# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile
