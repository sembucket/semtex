##############################################################################
# Remake linear algebra primitive library (liblinalp.a).
#
# Uses GNU make.
#
# RCS Information:
# $Id$
##############################################################################

#OPTIONS  = -O2 -mips2 -sopt
OPTIONS  = -g

LIB      = liblinalp.a

CFLAGS	 = $(OPTIONS)
ARFLAGS	 = cr
CPPFLAGS = -I.

#-----------------------------------------------------------------------------
# Create dependency list.
#-----------------------------------------------------------------------------

UTIL       = util

MEMORY     = memory

MATH       = xfill   xneg    xvneg   xsadd   xvadd   xssub   xvsub \
             xsmul   xvmul   xsdiv   xvrecp  xvdiv   xzero

OTHER      = xvabs   xvamax  xvpow   xvexp   xvlg10  xvlog   xvatan \
             xvatn2  xvcos   xvsin   xvsqrt  xrand   xfft

TRIAD      = xsvmvt  xsvpvt  xsvtsp  xsvtvm  xsvtvp  xsvvmt  xsvvpt \
	     xsvvtm  xsvvtp  xvvmvt  xvvpvt  xvvtvm  xvvtvp  xvvvtm

RELATIONAL = xseq    xsge    xsle    xslt    xsne

REDUCTION  = icount  ifirst  ixmax   ixmin   xsum

CONVERSION = vdble   xvfloa  vsngl   xbrev

MISC       = xscatr  xgathr  xramp   xcndst  xmask   xvpoly  xpolint \
	     xspline xmxm    xmxv    xmxva   fftdf 

ALL        = $(UTIL)       $(MEMORY)    $(MATH)       $(OTHER)   $(TRIAD) \
             $(RELATIONAL) $(REDUCTION) $(CONVERSION) $(MISC)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

#-----------------------------------------------------------------------------

$(LIB): $(LMD)

install:
	mv $(LIB) $(HOME)/lib; \
	cd $(HOME)/include; rm -f alplib.h; \
	ln -s $(HOME)/lib/linalp/alplib.h alplib.h

tar:
	tar cvf linalp.tar *

clean:
	rm *.o *.a *~ 
