##############################################################################
# Remake finite element routine library (libfem.a).
#
# Uses GNU make.
#
# $Id$
##############################################################################

ARCH = $(shell uname -s)

OPTIONS  = -g

FEMHOME  = ..
LIB      = libfem.a

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
endif

CC	 = cc
FC	 = f77
CFLAGS	 = $(OPTIONS)
FFLAGS   = $(OPTIONS)
ARFLAGS	 = cr
CPPFLAGS = -I$(FEMHOME)/include
RANLIB   =

# -------------- Compiler options for specific machines. ---------------------

ifeq ($(ARCH),IRIX)			# -- SGI IRIX 32-bit
CFLAGS     = $(DEFINES) -fullwarn -O2
FFLAGS     = $(DEFINES) -O2
endif

ifeq ($(ARCH),IRIX64)			# -- SGI IRIX 64-bit
CFLAGS     = $(DEFINES) -g -fullwarn
FFLAGS     = $(DEFINES) -g
endif

ifeq ($(ARCH),OSF1)			# -- DEC Alpha OSF1
FC	   = f90
CFLAGS     = $(DEFINES) -migrate -O4 -tune host
FFLAGS     = $(DEFINES) -migrate -O4 -tune host
endif

ifeq ($(ARCH),Linux)			# -- GNU/Linux system
CC         = gcc
FC         = g77
CFLAGS     = $(DEFINES) -O3
FFLAGS     = $(DEFINES) -O3
RANLIB     = ranlib
endif

# ----------------------------------------------------------------------------
# Create dependency list.
# ----------------------------------------------------------------------------

INITIAL    = initial

OPERATORS  = polyops   operators    norm     fftpack    fourier

MEMORY     = mapping   family

MESH       = RCM

ALL        = $(INITIAL) $(OPERATORS) $(MEMORY) $(MESH)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

# ----------------------------------------------------------------------------

# -- Make library, default action.

$(LIB): $(LMD)
ifeq ($(RANLIB),ranlib)
	$(RANLIB) $(LIB)
endif

# -- Install library & associated header files.

install:
	if test ! -d $(FEMHOME)/lib 		;		then	\
		mkdir $(FEMHOME)/lib 		; 		fi;	\
	if test ! -d $(FEMHOME)/lib/$(ARCH) 	;		then	\
		mkdir $(FEMHOME)/lib/$(ARCH)	; 		fi;	\
	if test -r $(LIB) 			;		then	\
		rm -f $(FEMHOME)/lib/$(ARCH)/$(LIB)	;		\
		cp $(LIB) $(FEMHOME)/lib/$(ARCH)	;	fi

	rm -f $(FEMHOME)/include/femdef.h;				\
	rm -f $(FEMHOME)/include/femlib.h;				\
	rm -f $(FEMHOME)/include/Femlib.h;				\

	cp femdef.h    $(FEMHOME)/include	;			\
	cp femlib.h    $(FEMHOME)/include	;			\
	cp Femlib.h    $(FEMHOME)/include	;			\

initial.o: defaults.h

testpoly: testpoly.c
	cc $(CFLAGS) -o testpoly testpoly.c -I$(FEMHOME)/include -I. \
	-L$(FEMHOME)/lib -lfem -lalp -lm -lmalloc

tar:
	tar cvf femlib.tar *

clean:
	rm *.o *.a *~ 
