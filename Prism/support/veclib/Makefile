#
# VECLIB - Library of Vector Routines
#

include $(PRISM)/system/defaults

LIB      = libvec.a
LIBDIR   = $(PRISM)/lib
INCDIR   = $(PRISM)/include
CPPFLAGS = -I$(PRISM)/support

VPATH    = math:triads:relational:reduction:conversion:misc:memory:blas:lapack

# These are the module names:

MATH    = xfill   xneg   xvneg  xvadd  xvsub  xsmul  xvmul  xvrecp  \
	  xvdiv   xzero  xsadd  xsdiv  xvabs  xvamax xvsin  xvcos   \
	  xvsqrt  xvlg10 xvexp  

TRIADS  = xsvtvm  xsvtvp  xsvvmt xsvvpt xsvvtm xsvvtp xvvtvm \
	  xvvtvp  xvvpvt  xvvmvt xsvtsp xvvvtm 

RELATE  = xseq    xsle    xslt   xsne

REDUCE  = xsum    ixmax   ixmin  icount ifirst xvpoly

CONVERT = xvfloa  vdble   vsngl  xvrand xvreal xvimag xvcmplx xbrev

MISC    = xscatr  xgathr  xramp  xcndst xmask  polint xfft   \
          mxm     mxv     mxva   fftdf  spline veclib 

ifneq ($(ARCH),i860)
override MISC := $(MISC) dclock	
endif

MEMORY  = xvector xmatrix 

ifndef LAPACK
LAPACK  = dpptrf  dpotf2  dpptrs dppcon \
	  spptrf  spotf2  spptrs sppcon \
	  dpbtrf  dpbtf2  dpbtrs dpbcon \
	  spbtrf  spbtf2  spbtrs spbcon \
	  dgetrf  dgetf2  dgetrs dgecon dgesv \
	  dgeev   dgebak  dgebal dgehrd \
	  dlaswp  dlamch  dlamch dlacon dlatbs dlatps dlabad dlange \
	  dlascl  dlacpy  dorghr dhseqr dtrevc dlapy2 dlartg dgehd2 \
	  dlahrd  dlarfb  dlassq dlaset dlahqr dlanhs dlarfx dlarfg \
	  dlaln2  dorgqr  dlarf  dlanv2 dladiv dorg2r dlarft \
	  slaswp  slamch  slamch slacon slatbs  slatps  slabad \
	  drscl   srscl   ilaenv  
else
override LAPACK :=
endif

ifndef BLAS
BLAS    = dasum   daxpy   dcopy   ddot   dgemm  dgemv  dger   dnrm2  \
	  dscal   dspr    dswap   dsyr   dsyrk  dtbsv  dtpsv  dtrsm  \
	  idamax  lsame   xerbla
else
override BLAS := 
endif

ALL     = $(MATH) $(TRIADS) $(RELATE) $(REDUCE) $(CONVERT) $(MISC) \
          $(MEMORY) $(LAPACK) $(BLAS)
HEADER  = veclib.h blas.h lapack.h

# Generate rules for building the modules

LMOD    = $(foreach routine, $(ALL), $(LIB)($(routine).o))
OBJS    = $(foreach routine, $(ALL), $(routine).o)
TARFILE = $(patsubst %.a, %.tar, $(LIB))

# -----------------------------------------------------------

$(LIB): $(LMOD)
	$(RANLIB) $(LIB)

install: $(LIB)
	cp $(LIB) $(LIBDIR)/$(LIB)
	if [ ! -d $(INCDIR)/veclib ]; \
	then mkdir $(INCDIR)/veclib; \
	else rm -f $(INCDIR)/veclib/*; \
	fi
	cp $(HEADER) $(INCDIR)/veclib/.

uninstall:
	-rm -f  $(LIBDIR)/$(LIB)
	-rm -rf $(INCDIR)/veclib


tar:
	-rm -f $(TARFILE)
	tar cvf $(TARFILE) README Makefile veclib.h blas.h lapack.h \
		blas conversion lapack math memory misc/*.c misc/*.f \
		reduction relational triads

tidy:
	-rm -f *~ \#*\#
clean: tidy
	-rm -f *.o
empty: clean
	-rm -f $(LIB) $(TARFILE)

# The Cray compiler with the given options generates bugs in the
# following files:

ifeq ($(ARCH),CRAY)
spline.o: spline.c
	$(CC) -c $<
endif
