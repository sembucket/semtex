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
CFLAGS	 = $(OPTIONS)
ARFLAGS	 = cr
RANLIB   =

# ----------- Compiler options for specific machines: ------------------------

ifeq ($(ARCH),IRIX)			# -- SGI IRIX 5.3
CPPFLAGS   = -I../include
CFLAGS     = -fullwarn -O2 -sopt -mips2
endif

ifeq ($(ARCH),IRIX64)			# -- SGI IRIX 64-bit
CPPFLAGS   = -I../include
CFLAGS     = -fullwarn -mips4 -O3
endif

ifeq ($(ARCH),OSF1)			# -- DEC Alpha OSF1
CPPFLAGS   = -I../include
CC         = cc
CFLAGS     = -migrate -fast -O4 -tune host
#CC	   = kcc
#CFLAGS    = -ckapargs='-optimize=5 -roundoff=3 -scaleropt=3'
endif

ifeq ($(ARCH),Linux)			# -- GNU/Linux system
CPPFLAGS   = -I../include
CC         = gcc
CFLAGS     = -O3
RANLIB     = ranlib
endif

ifeq ($(ARCH),UNIX_System_V)            # -- Fujitsu VPP300
CPPFLAGS   = -I../include
CC         = vcc
#CFLAGS     = -Wv,-ad
CFLAGS     = -Wv,-m3,-Om
RANLIB     = 
endif

ifeq ($(ARCH),SUPER-UX)                 # -- NEC SX-4
CPPFLAGS   = -I../include
CC         = cc
CFLAGS     = -pvctl,fullmsg -c -O nooverlap
RANLIB     = 
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
	if test ! -d    $(FEMHOME)/lib          ;               then    \
		mkdir   $(FEMHOME)/lib          ;               fi;     \
	if test ! -d    $(FEMHOME)/lib/$(ARCH)  ;               then    \
		mkdir   $(FEMHOME)/lib/$(ARCH)  ;               fi;     \
	if test -r      $(LIB)                  ;               then    \
		rm -f   $(FEMHOME)/lib/$(ARCH)/$(LIB)   ;               \
		cp $(LIB) $(FEMHOME)/lib/$(ARCH)        ;       fi

	if test ! -d   $(FEMHOME)/include       ;               then    \
		mkdir $(FEMHOME)/include       ;               fi

	rm -f $(FEMHOME)/include/alplib.h;      \
	rm -f $(FEMHOME)/include/Veclib.h;      \
	rm -f $(FEMHOME)/include/Utility.h;     \
	rm -f $(FEMHOME)/include/Blas.h;        \
	rm -f $(FEMHOME)/include/Lapack.h

	cp alplib.h    $(FEMHOME)/include       ;                       \
	cp Veclib.h    $(FEMHOME)/include       ;                       \
	cp Utility.h   $(FEMHOME)/include       ;                       \
	cp Blas.h      $(FEMHOME)/include       ;                       \
	cp Lapack.h    $(FEMHOME)/include

tar:
	tar cvf alplib.tar *

clean:
	rm *.o *.a *~ 

ifeq ($(ARCH),IRIX) # -- special compile to avoid bug
xbrev.o: xbrev.c
	cc $(CPPFLAGS) -c -O1 -mips2 xbrev.c
endif
