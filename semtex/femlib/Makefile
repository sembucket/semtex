##############################################################################
# Remake finite element routine library (libfem.a).  Gnu make required.
#
# $Id$
##############################################################################

SEM    = ..
SEMSRC = $(SEM)/src
ARCH   = $(shell uname -s)

ifdef MPI
  LIB = libfem_mp.a
else
  LIB = libfem.a
endif

ifdef DEBUG_FFT
  DEFINES += -DDEBUG_FFT
endif

# ----------------------------------------------------------------------------
# Create dependency list.
#

INITIAL    = initial
#INITIAL    = manager tree

OPERATORS  = polyops	operators	matops     polylib     filter

FFT	   = temfftd	fftpack		fourier

MEMORY     = mapping	family          message

MESH       = RCM

IMPORT     = recipes zbesj mach xerbla

ALL        = $(INITIAL) $(OPERATORS) $(FFT) $(MEMORY) $(MESH) $(IMPORT)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

# ----------------------------------------------------------------------------
# Make library, default action.
#
$(LIB): $(LMD)
#ifdef $(RANLIB) # -- Don't know why this won't work yet...
#$(RANLIB) $(LIB)
#endif
ifeq ($(ARCH),Darwin) # -- So use workaround.
	ranlib $(LIB)
endif
ifeq ($(ARCH),Linux-alpha)
	ranlib $(LIB)
endif
ifeq ($(ARCH),Linux-iX86)
	ranlib $(LIB)
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
		rm -f $(SEM)/lib/$(ARCH)/$(LIB)	;			\
		cp $(LIB) $(SEM)/lib/$(ARCH)	;		fi

	rm -f $(SEM)/include/cfemdef	;			\
	rm -f $(SEM)/include/cfemlib    ;			\
	rm -f $(SEM)/include/cnr77      ;			\
	rm -f $(SEM)/include/cpolylib   ;			\
	rm -f $(SEM)/include/femlib_h	;			

	cp cfemdef	$(SEM)/include	;			\
	cp cfemlib	$(SEM)/include	;			\
	cp cnr77	$(SEM)/include	;			\
	cp cpolylib	$(SEM)/include  ;			\
	cp femlib_h	$(SEM)/include

# ----------------------------------------------------------------------------
# Ensure any new default parser variables are picked up.
#
initial.o: cdefaults

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
