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

SEMHOME = ..
LIB     = libalp.a

# ----------------------------------------------------------------------------
# Generic settings; see specific settings for specific machines below.
OPTIONS  = -g

CC       = cc
OPTIONS  = -g
CFLAGS	 = $(OPTIONS)
CPPFLAGS = -I$(SEMHOME)/include
ARFLAGS	 = cr
RANLIB   =

# ----------- Compiler options for specific machines: ------------------------

ifeq ($(ARCH),IRIX)			# -- SGI IRIX 32-bit.
CFLAGS     = -fullwarn -O2 -sopt
endif

ifeq ($(ARCH),IRIX64)			# -- SGI IRIX 64-bit
CPPFLAGS   = -I../include
CFLAGS     = -fullwarn -mips4 -O3
endif

ifeq ($(ARCH),OSF1)			# -- DEC Alpha OSF1
CFLAGS     = -migrate -fast -O4 -tune host
#CC	   = kcc
#CFLAGS    = -ckapargs='-optimize=5 -roundoff=3 -scaleropt=3'
endif

ifeq ($(ARCH),Linux)			# -- GNU/Linux system
CC         = gcc
CFLAGS     = -O3
RANLIB     = ranlib
endif

ifeq ($(ARCH),UNIX_System_V)            # -- Fujitsu VPP300
CC         = vcc
#CFLAGS     = -Wv,-ad
CFLAGS     = -Wv,-m3,-Om
endif

ifeq ($(ARCH),SUPER-UX)                 # -- NEC SX-4
CFLAGS     = -pvctl,fullmsg -O nooverlap -hfloat0
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

MISC       = xscatr      xgathr      xgathr_scatr	\
	     xscatr_sum  xgathr_sum  xgathr_scatr_sum	\
	     xramp       xcndst      xmask       	\
	     xvpoly      xpolint    			\
	     xspline     xsplquad			\
	     xmxm        xmxv        xmxva	  	\
	     xvvtvvtp    xsvvttvp

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
	if test ! -d    $(SEMHOME)/lib          ;               then    \
		mkdir   $(SEMHOME)/lib          ;               fi;     \
	if test ! -d    $(SEMHOME)/lib/$(ARCH)  ;               then    \
		mkdir   $(SEMHOME)/lib/$(ARCH)  ;               fi;     \
	if test -r      $(LIB)                  ;               then    \
		rm -f   $(SEMHOME)/lib/$(ARCH)/$(LIB)   ;               \
		cp $(LIB) $(SEMHOME)/lib/$(ARCH)        ;       fi

	if test ! -d   $(SEMHOME)/include       ;               then    \
		mkdir $(SEMHOME)/include       ;               fi

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
	rm *.o *.a *~ 
