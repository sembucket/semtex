##############################################################################
# Remake finite element routine library (libfem.a).
#
# Uses GNU make.
#
# RCS Information:
# $Id$
##############################################################################

ARCH = $(shell uname -s)

OPTIONS  = -O2 -sopt

FEMHOME  = ../..
LIBHOME  = $(shell pwd)
LIB      = libfem.a

CFLAGS	 = $(OPTIONS)
FFLAGS   = $(OPTIONS)
ARFLAGS	 = cr
CPPFLAGS = -I../../include

# -------------- compiler options for specific machines. ---------------------

ifeq ($(ARCH),IRIX)
CFLAGS     = -mips2 -fullwarn $(OPTIONS)
FFLAGS     = -mips2  $(OPTIONS)
endif

ifeq ($(ARCH),IRIX64)
CFLAGS     = -fullwarn $(OPTIONS)
FFLAGS     = $(OPTIONS)
endif

#-----------------------------------------------------------------------------
# Create dependency list.
#-----------------------------------------------------------------------------

INITIAL    = initial

OPERATORS  = polyops   operators    norm

MESH       = RCM

ALL        = $(INITIAL) $(OPERATORS) $(MESH)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

#-----------------------------------------------------------------------------

$(LIB): $(LMD)

install:
	cd $(FEMHOME)/lib; rm -f $(LIB);			\
	ln -s $(LIBHOME)/$(LIB) $(LIB);				\
	cd $(HOME)/lib; rm -f $(LIB);				\
	ln -s $(LIBHOME)/$(LIB) $(LIB);				\
	cd $(HOME)/include; 					\
	rm -f femdef.h femlib.h Femlib.h;                     	\
	ln -s $(LIBHOME)/femdef.h    femdef.h;			\
	ln -s $(LIBHOME)/femlib.h    femlib.h;			\
	ln -s $(LIBHOME)/Femlib.h    Femlib.h;

	cd $(FEMHOME)/include; 					\
	rm -f femdef.h femlib.h Femlib.h;			\
	ln -s $(LIBHOME)/femdef.h    femdef.h;			\
	ln -s $(LIBHOME)/femlib.h    femlib.h;			\
	ln -s $(LIBHOME)/Femlib.h    Femlib.h

initial.y: defaults.h

testpoly: testpoly.c
	cc $(CFLAGS) -o testpoly testpoly.c -I$(FEMHOME)/include -I. \
	-L$(FEMHOME)/lib -lfem -lalp -lm -lmalloc

tar:
	tar cvf femlib.tar *

clean:
	rm *.o *.a *~ 
