##############################################################################
# Remake linear algebra primitive library (libalp.a), and maintain header
# files for Veclib, BLAS, & LAPACK calls.
#
# Uses GNU make.
#
# RCS Information:
# $Id$
##############################################################################

ARCH    = $(shell uname -s)

FEMHOME = ..
LIB     = libalp.a

# ----------------------------------------------------------------------------
# Generic settings; see specific settings for specific machines below.
OPTIONS  = -g

CC       = cc
FC	 = f77
CFLAGS	 = $(OPTIONS)
FFLAGS   = $(OPTIONS)
ARFLAGS	 = cr
RANLIB   =

# ----------- Compiler options for specific machines: ------------------------

ifeq ($(ARCH),IRIX)			# -- SGI IRIX 5.3
CFLAGS     = -fullwarn -O2 -sopt
FFLAGS     = -O2 -sopt
endif

ifeq ($(ARCH),IRIX64)			# -- SGI IRIX 64-bit
CFLAGS     = -fullwarn -g
FFLAGS     = -g
endif

ifeq ($(ARCH),OSF1)			# -- DEC Alpha OSF1
FC	   = f90
CFLAGS     = -migrate -O4 -tune host
FFLAGS     = -migrate -O4 -tune host
endif

ifeq ($(ARCH),Linux)			# -- GNU/Linux system
CC         = gcc
FC         = g77
CFLAGS     = -O3
FFLAGS     = -O3
RANLIB     = ranlib
endif

# ----------------------------------------------------------------------------
# Create dependency list.
# ----------------------------------------------------------------------------

UTIL       = util

MEMORY     = memory

MATH       = xcopy   xfill   xneg    xvneg   xsadd   xvadd   xssub \
             xvsub   xsmul   xvmul   xsdiv   xvrecp  xvdiv   xzero   

OTHER      = xvabs   xvamax  xvpow   xvexp   xvlg10  xvlog   xvatan \
             xvatn2  xvcos   xvsin   xvsqrt  xrand   xfft    xvhypot

TRIAD      = xsvmvt  xsvpvt  xsvtsp  xsvtvm  xsvtvp  xsvvmt  xsvvpt \
	     xsvvtm  xsvvtp  xvvmvt  xvvpvt  xvvtvm  xvvtvp  xvvvtm

RELATIONAL = xseq    xsge    xsle    xslt    xsne

REDUCTION  = icount  ifirst  ixmax   ixmin   xsum    lany    lxsame

CONVERSION = vdble   xvfloa  vsngl   xbrev

MISC       = xscatr          xgathr          xgathr_scatr	\
	     xscatr_sum      xgathr_sum      xgathr_scatr_sum	\
	     xramp   xcndst  xmask   xvpoly  xpolint		\
	     xspline xmxm    xmxv    xmxva   fftdf

ALL        = $(UTIL)       $(MEMORY)    $(MATH)       $(OTHER)   $(TRIAD) \
             $(RELATIONAL) $(REDUCTION) $(CONVERSION) $(MISC)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

# ----------------------------------------------------------------------------

# -- Make library, default action.

$(LIB): $(LMD)
ifeq ($(RANLIB),ranlib)
	$(RANLIB) $(LIB)
endif

# -- Install library & associated header files.

install:
	if test ! -d	$(HOME)/lib/$(ARCH) 	;		then	\
		mkdir	$(HOME)/lib/$(ARCH)	;		fi;	\
	if test ! -d	$(FEMHOME)/lib		;		then	\
		mkdir	$(FEMHOME)/lib		; 		fi;	\
	if test ! -d	$(FEMHOME)/lib/$(ARCH) 	;		then	\
		mkdir	$(FEMHOME)/lib/$(ARCH)	; 		fi;	\
	if test -r	$(LIB) 			;		then	\
		rm -f	$(FEMHOME)/lib/$(ARCH)/$(LIB)	;		\
		cp $(LIB) $(FEMHOME)/lib/$(ARCH)	;		\
		rm -f	$(HOME)/lib/$(ARCH)/$(LIB)	;		\
		cp $(LIB) $(HOME)/lib/$(ARCH)		;	fi

	if test ! -d   $(HOME)/include		;		then	\
		 mkdir $(HOME)/include		;		fi;	\
	if test ! -d   $(FEMHOME)/include	;		then	\
		 mkdir $(FEMHOME)/include	;		fi

	rm -f $(HOME)/include/alplib.h	$(FEMHOME)/include/alplib.h;	\
	rm -f $(HOME)/include/Veclib.h	$(FEMHOME)/include/Veclib.h;	\
	rm -f $(HOME)/include/Utility.h $(FEMHOME)/include/Utility.h;	\
	rm -f $(HOME)/include/Blas.h 	$(FEMHOME)/include/Blas.h;	\
	rm -f $(HOME)/include/Lapack.h	$(FEMHOME)/include/Lapack.h

	cp alplib.h    $(HOME)/include		;			\
	cp alplib.h    $(FEMHOME)/include	;			\
	cp Veclib.h    $(HOME)/include		;			\
	cp Veclib.h    $(FEMHOME)/include	;			\
	cp Utility.h   $(HOME)/include		;			\
	cp Utility.h   $(FEMHOME)/include	;			\
	cp Blas.h      $(HOME)/include		;			\
	cp Blas.h      $(FEMHOME)/include	;			\
	cp Lapack.h    $(HOME)/include		;			\
	cp Lapack.h    $(FEMHOME)/include

tar:
	tar cvf alplib.tar *

clean:
	rm *.o *.a *~ 

ifeq ($(ARCH),IRIX) # -- special compile to avoid bug
xbrev.o: xbrev.c
	cc -c -O1 -mips2 xbrev.c
endif
