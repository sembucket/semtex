##############################################################################
# Remake linear algebra primitive library (libalp.a), and maintain header
# files for Veclib, BLAS, & LAPACK calls.
#
# $Id$
##############################################################################

SEM    = ..
SEMSRC = $(SEM)/src
LIB    = libalp.a

# ----------------------------------------------------------------------------
# Create dependency list.
# ----------------------------------------------------------------------------

UTIL       = util

MEMORY     = memory

MATH       = xcopy   xfill   xneg    xvneg   xsadd   xvadd   xssub  \
             xvsub   xsmul   xvmul   xsdiv   xvrecp  xvdiv   xzero  \
	     xspow   xvsgn

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
	     xvpoly      xpolint    	xclip			\
	     xspline     xsplquad	xiclip			\
	     xmxv        xmxva	  				\
	     xvvtvvtp    xsvvttvp	 xvvtvvtm

FORTRAN    = xmxm							

ALL        = $(UTIL)       $(MEMORY)    $(MATH)       $(OTHER)   $(TRIAD) \
             $(RELATIONAL) $(REDUCTION) $(CONVERSION) $(MISC)    $(FORTRAN)    
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
	if test ! -d    $(SEM)/lib          ;               then    \
		mkdir   $(SEM)/lib          ;               fi;     \
	if test ! -d    $(SEM)/lib/$(ARCH)  ;               then    \
		mkdir   $(SEM)/lib/$(ARCH)  ;               fi;     \
	if test -r      $(LIB)                  ;               then    \
		rm -f   $(SEM)/lib/$(ARCH)/$(LIB)   ;               \
		cp $(LIB) $(SEM)/lib/$(ARCH)        ;       fi

	if test ! -d   $(SEM)/include       ;               then    \
		mkdir $(SEM)/include        ;               fi

	rm -f $(SEM)/include/alplib.h;      \
	rm -f $(SEM)/include/Veclib.h;      \
	rm -f $(SEM)/include/Utility.h;     \
	rm -f $(SEM)/include/Blas.h;        \
	rm -f $(SEM)/include/Lapack.h

	cp alplib.h    $(SEM)/include       ;                       \
	cp Veclib.h    $(SEM)/include       ;                       \
	cp Utility.h   $(SEM)/include       ;                       \
	cp Blas.h      $(SEM)/include       ;                       \
	cp Lapack.h    $(SEM)/include


