##############################################################################
# Remake finite element routine library (libfem.a).
#
# Uses GNU make.
#
# RCS Information:
# $Id$
##############################################################################

ARCH = $(shell uname -s)

OPTIONS  = -O2
#OPTIONS  = -g

LIB      = libfem.a

CFLAGS	 = $(OPTIONS)
FFLAGS   = $(OPTIONS)
ARFLAGS	 = cr
CPPFLAGS = -I../../include

# ----------- compiler options for specific machines.

ifeq ($(ARCH),IRIX)
CFLAGS     = -mips2 -sopt -fullwarn $(OPTIONS)
FFLAGS     = -mips2 -sopt $(OPTIONS)
endif

ifeq ($(ARCH),IRIX64)
CFLAGS     = -sopt -fullwarn $(OPTIONS)
FFLAGS     = -sopt $(OPTIONS)
endif


#-----------------------------------------------------------------------------
# Create dependency list.
#-----------------------------------------------------------------------------

INITIAL    = initial

OPERATORS  = polyops   operators

MESH       = RCM

ALL        = $(INITIAL) $(OPERATORS) $(MESH)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

#-----------------------------------------------------------------------------

$(LIB): $(LMD)

install:
	mv $(LIB) ..; \
	cd ../../include; rm -f femlib.h; ln -s ../lib/femlib/femlib.h femlib.h

initial.y: defaults.h

testpoly: testpoly.c
	cc $(CFLAGS) -o testpoly testpoly.c -I$(HOME)/include -I. \
	-L$(HOME)/lib -L. -lfem -llinalp -lm -lmalloc

tar:
	tar cvf femlib.tar *

clean:
	rm *.o *.a *~ 
