##############################################################################
# Remake finite element routine library (libfem.a).  Gnu make required.
#
# $Id$
##############################################################################

SEM = ..
SEMSRC = $(SEM)/src

ifdef MPI
  LIB = libfem_mp.a
else
  LIB = libfem.a
endif

# ----------------------------------------------------------------------------
# Create dependency list.
#

#INITIAL    = initial
INITIAL    = manager tree

OPERATORS  = polyops	operators	matops          filter

FFT	   = temfftd	fftpack		fourier

MEMORY     = mapping	family          message

MESH       = RCM

IMPORT     = recipes zbesj mach

ALL        = $(INITIAL) $(OPERATORS) $(FFT) $(MEMORY) $(MESH) $(IMPORT)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

# ----------------------------------------------------------------------------
# Make library, default action.
#
$(LIB): $(LMD)
ifeq ($(RANLIB),ranlib)
	$(RANLIB) $(LIB)
endif

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(SEMSRC)/Makefile

# ----------------------------------------------------------------------------
# Install library & associated header files.
#
install:
	if test ! -d $(SEM)/lib 		;		then	\
		mkdir $(SEM)/lib 		; 		fi;	\
	if test ! -d $(SEM)/lib/$(ARCH) 	;		then	\
		mkdir $(SEM)/lib/$(ARCH)	; 		fi;	\
	if test -r $(LIB) 			;		then	\
		rm -f $(SEM)/lib/$(ARCH)/$(LIB)	;		\
		cp $(LIB) $(SEM)/lib/$(ARCH)	;	fi

	rm -f $(SEM)/include/femdef.h;				\
	rm -f $(SEM)/include/femlib.h;				\
	rm -f $(SEM)/include/Femlib.h;				\
	rm -f $(SEM)/include/nr77.h;				\
	rm -r $(SEM)/include/polylib.h

	cp femdef.h    $(SEM)/include	;			\
	cp femlib.h    $(SEM)/include	;			\
	cp Femlib.h    $(SEM)/include	;			\
	cp nr77.h      $(SEM)/include	;			\
	cp polylib.h   $(SEM)/include

# ----------------------------------------------------------------------------
# Ensure any new default parser variables are picked up.
#
initial.o: defaults.h

# ----------------------------------------------------------------------------
# Special machine-dependent compilations.
#
ifeq ($(ARCH),OSF1)
temfftd.o: temfftd.F
	f90 -c -fast -inline speed -tune host temfftd.F -O5 -unroll 2
temffts.o: temffts.F
	f90 -c -fast -inline speed -tune host temffts.F -O5 -unroll 2
endif

# ----------------------------------------------------------------------------
# Rules to help with .F routines.
#
matops.o:  matops.F
temfftd.o: temfftd.F
