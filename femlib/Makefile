##############################################################################
# Remake finite element routine library (libfem.a).
#
# Uses GNU make.
#
# RCS Information:
# $Id$
##############################################################################

OPTIONS  = -g

LIB      = libfem.a

CFLAGS	 = $(OPTIONS)
ARFLAGS	 = cr
CPPFLAGS = -I$(HOME)/include -I../../include

#-----------------------------------------------------------------------------
# Create dependency list.
#-----------------------------------------------------------------------------

INITIAL    = initial

ALL        = $(INITIAL)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

#-----------------------------------------------------------------------------

$(LIB): $(LMD)

install:
	mv $(LIB) ..; \

initial.y: defaults.h

tar:
	tar cvf femlib.tar *

clean:
	rm *.o *.a *~ 
