##############################################################################
# Remake linear algebra primitive library (libalp.a), and maintain header
# files for Veclib, BLAS, & LAPACK calls.
#
# Uses GNU make.
#
# RCS Information:
# $Id$
##############################################################################

SEMHOME = ..
LIB     = libalp.a

ARCH = $(shell uname -s)
MACH = $(shell uname -m)
ifeq ($(ARCH),Linux)
  ifeq ($(MACH),alpha)
    ARCH = Linux-alpha
  endif
  ifeq ($(MACH),i586)
    ARCH = Linux-i586
  endif
endif

# ----------------------------------------------------------------------------
# Generic settings; see specific settings for specific machines below.

OPTIONS  = -g
CC       = cc
FC	 = f77
OPTIONS  = -g
CFLAGS	 = $(OPTIONS)
CPPFLAGS = -I$(SEMHOME)/include
ARFLAGS	 = cr
RANLIB   =

# ----------- Compiler options for specific machines: ------------------------

ifeq ($(ARCH),IRIX)				# -- SGI IRIX 32-bit.
CFLAGS = -fullwarn -O2
FFLAGS = -O2
endif

ifeq ($(ARCH),IRIX64)				# -- SGI IRIX 64-bit
CFLAGS = -mips4 -64 -r10000 -Ofast=IP30 -OPT:roundoff=3
FFLAGS = -mips4 -64 -r10000 -Ofast=IP30 -OPT:roundoff=3
endif


ifeq ($(ARCH),OSF1)				# -- DEC Alpha OSF1
FC = f90
CC = cc
CFLAGS = -fast -inline speed -fp_reorder -speculate all -tune host -O4
FFLAGS = -fast -tune host -O5
endif

ifeq ($(ARCH),Linux-i586)			# -- Linux Intel system
CC = gcc
FC = g77
CFLAGS = -O3
CFLAGS = -O3
RANLIB = ranlib
endif

ifeq ($(ARCH),Linux-alpha)			# -- Linux Alpha system
CC = ccc
FC = fort
CFLAGS = -arch host -tune host -fast -O4
FFLAGS = -arch host -tune host -fast -O4
RANLIB = ranlib
endif

ifeq ($(ARCH),UNIX_System_V)			# -- Fujitsu VPP300
CC = vcc
FC = frt
CFLAGS = -Wv,-m3,-Om -Kpopt
FFLAGS = -Wv,-m3,-te
endif

ifeq ($(ARCH),SUPER-UX)                 	# -- NEC SX-4
FC       = f90
DEFINES += -D_BUILTIN_
CFLAGS   = $(DEFINES) -pvctl,fullmsg -Onooverlap -hfloat0 -hnostkchk
FFLAGS   = $(DEFINES) -Wf'-pvctl,fullmsg -ptr byte -b' -float0 -ew
endif

# ----------------------------------------------------------------------------
# Create dependency list.
# ----------------------------------------------------------------------------

UTIL       = util

MEMORY     = memory

MATH       = xcopy   xfill   xneg    xvneg   xsadd   xvadd   xssub  \
             xvsub   xsmul   xvmul   xsdiv   xvrecp  xvdiv   xzero  \
	     xspow

OTHER      = xvabs   xvamax  xvpow   xvexp   xvlg10  xvlog   xvatan \
             xvatn2  xvcos   xvsin   xvsqrt  xrand   xvhypot

TRIAD      = xsvmvt  xsvpvt  xsvtsp  xsvtvm  xsvtvp  xsvvmt  xsvvpt \
	     xsvvtm  xsvvtp  xvvmvt  xvvpvt  xvvtvm  xvvtvp  xvvvtm \
	     xsvvtt

RELATIONAL = xseq    xsge    xsle    xslt    xsne

REDUCTION  = icount  ifirst  ixmax   ixmin   xsum    lany    lxsame

CONVERSION = vdble   xvfloa  vsngl   xbrev

MISC       = xscatr      xgathr      	xgathr_scatr		\
	     xscatr_sum  xgathr_sum  	xgathr_scatr_sum	\
	     xramp       xcndst      	xmask       		\
	     xvpoly      xpolint    				\
	     xspline     xsplquad				\
	     xmxv        xmxva	  				\
	     xvvtvvtp    xsvvttvp	 xvvtvvtm

FORTRAN    = xmxm							

ALL        = $(UTIL)       $(MEMORY)    $(MATH)       $(OTHER)   $(TRIAD) \
             $(RELATIONAL) $(REDUCTION) $(CONVERSION) $(MISC)    $(FORTRAN)    
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

# ----------------------------------------------------------------------------

# -- Make library, default action.

$(LIB): $(LMD)
ifeq ($(RANLIB),ranlib)
	$(RANLIB) $(LIB)
endif

# -- Install library & associated header files.

install:
	if test ! -d    $(SEMHOME)/lib          ;               then    \
		mkdir   $(SEMHOME)/lib          ;               fi;     \
	if test ! -d    $(SEMHOME)/lib/$(ARCH)  ;               then    \
		mkdir   $(SEMHOME)/lib/$(ARCH)  ;               fi;     \
	if test -r      $(LIB)                  ;               then    \
		rm -f   $(SEMHOME)/lib/$(ARCH)/$(LIB)   ;               \
		cp $(LIB) $(SEMHOME)/lib/$(ARCH)        ;       fi

	if test ! -d   $(SEMHOME)/include       ;               then    \
		mkdir $(SEMHOME)/include        ;               fi

	rm -f $(SEMHOME)/include/alplib.h;      \
	rm -f $(SEMHOME)/include/Veclib.h;      \
	rm -f $(SEMHOME)/include/Utility.h;     \
	rm -f $(SEMHOME)/include/Blas.h;        \
	rm -f $(SEMHOME)/include/Lapack.h

	cp alplib.h    $(SEMHOME)/include       ;                       \
	cp Veclib.h    $(SEMHOME)/include       ;                       \
	cp Utility.h   $(SEMHOME)/include       ;                       \
	cp Blas.h      $(SEMHOME)/include       ;                       \
	cp Lapack.h    $(SEMHOME)/include

tar:
	tar cvf alplib.tar *

clean:
	rm *.o *.a *.O *~ 
