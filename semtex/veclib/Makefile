##############################################################################
# Remake linear algebra primitive library (libalp.a).
#
# Uses GNU make.
#
# RCS Information:
# $Id$
##############################################################################

ARCH    = $(shell uname -s)

FEMHOME = ../..
LIBHOME = $(shell pwd)
LIB     = libalp.a

OPTIONS  = -O2

CC       = cc
CFLAGS	 = $(OPTIONS)
ARFLAGS	 = cr

# ----------- Compiler options for specific machines: ------------------------

ifeq ($(ARCH),IRIX)
CFLAGS     = -mips2 -sopt -fullwarn $(OPTIONS)
endif

ifeq ($(ARCH),IRIX64)
CFLAGS     = -sopt -fullwarn $(OPTIONS)
endif

#-----------------------------------------------------------------------------
# Create dependency list.
#-----------------------------------------------------------------------------

UTIL       = util

MEMORY     = memory

MATH       = xcopy   xfill   xneg    xvneg   xsadd   xvadd   xssub \
             xvsub   xsmul   xvmul   xsdiv   xvrecp  xvdiv   xzero   

OTHER      = xvabs   xvamax  xvpow   xvexp   xvlg10  xvlog   xvatan \
             xvatn2  xvcos   xvsin   xvsqrt  xrand   xfft    xvhypot

TRIAD      = xsvmvt  xsvpvt  xsvtsp  xsvtvm  xsvtvp  xsvvmt  xsvvpt \
	     xsvvtm  xsvvtp  xvvmvt  xvvpvt  xvvtvm  xvvtvp  xvvvtm

RELATIONAL = xseq    xsge    xsle    xslt    xsne

REDUCTION  = icount  ifirst  iany    ixmax   ixmin   xsum

CONVERSION = vdble   xvfloa  vsngl   xbrev

MISC       = xscatr  xgathr  xramp   xcndst  xmask   xvpoly  xpolint \
	     xspline xmxm    xmxv    xmxva   fftdf 

ALL        = $(UTIL)       $(MEMORY)    $(MATH)       $(OTHER)   $(TRIAD) \
             $(RELATIONAL) $(REDUCTION) $(CONVERSION) $(MISC)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

#-----------------------------------------------------------------------------

$(LIB): $(LMD)

install:
	cd $(FEMHOME)/lib; rm -f $(LIB);			\
	ln -s $(LIBHOME)/$(LIB) $(LIB)				\
	cd $(HOME)/lib; rm -f $(LIB);				\
	ln -s $(LIBHOME)/$(LIB) $(LIB);				\
	cd $(HOME)/include; rm -f alplib.h Veclib.h Util.h;	\
	ln -s $(LIBHOME)/alplib.h alplib.h;  			\
	ln -s $(LIBHOME)/Veclib.h Veclib.h;			\
	ln -s $(LIBHOME)/Utility.h Utility.h;				
	cd $(FEMHOME)/include; rm -f alplib.h Veclib.h Util.h;	\
	ln -s $(LIBHOME)/alplib.h alplib.h;  			\
	ln -s $(LIBHOME)/Veclib.h Veclib.h;			\
	ln -s $(LIBHOME)/Utility.h Utility.h;

tar:
	tar cvf alplib.tar *

clean:
	rm *.o *.a *~ 
