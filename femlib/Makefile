##############################################################################
# Remake finite element routine library (libfem.a).  Gnu make required.
#
# $Id$
##############################################################################

ARCH = $(shell uname -s)

SEMHOME  = ..
LIB      = libfem.a

DEFINES  =
ifdef DEBUG
  DEFINES += -DDEBUG
endif

CC	 = cc
CXX	 = CC
FC	 = f77
OPTIONS  = -g
CFLAGS	 = $(OPTIONS)
FFLAGS   = $(OPTIONS)
ARFLAGS	 = cr
CPPFLAGS = -I$(SEMHOME)/include
RANLIB   =

# -------------- Compiler options for specific machines. ---------------------

ifeq ($(ARCH),IRIX)			# -- SGI IRIX 32-bit
CFLAGS     = $(DEFINES) -fullwarn -O2
FFLAGS     = $(DEFINES) -O2
endif

ifeq ($(ARCH),IRIX64)			# -- SGI IRIX 64-bit
CFLAGS     = $(DEFINES) -fullwarn -mips4 -O2
FFLAGS     = $(DEFINES) -mips4 -O2
endif

ifeq ($(ARCH),OSF1)			# -- DEC Alpha OSF1
FC	   = f90
CXX	   = cxx
CC         = cc
CFLAGS     = $(DEFINES) -fast -inline speed \
              -tune host -O4
FFLAGS     = $(DEFINES) -fast -inline speed -tune host -O5 -unroll 2
ifdef MPI
  DEFINES  += -DMPI
  CPPFLAGS += -I/usr/local/mpi/include
  CFLAGS   += -pthread
  FFLAGS   += -D_THREAD_SAFE
endif
endif

ifeq ($(ARCH),Linux)			# -- GNU/Linux system
CC         = gcc
FC         = g77
CXX	   = g++
CFLAGS     = $(DEFINES) -O3
FFLAGS     = $(DEFINES) -O3
RANLIB     = ranlib
endif

ifeq ($(ARCH),UNIX_System_V)            # -- Fujitsu VPP300
CC         = vcc
FC         = frt
CFLAGS     = -Wv,-m3,-Om,-te
FFLAGS     = -Wv,-m3,-te
endif

ifeq ($(ARCH),SUPER-UX)                 # -- NEC SX-4
CFLAGS     = $(DEFINES) -pvctl,fullmsg -hfloat0
FFLAGS     = $(DEFINES) -Wf'-pvctl,fullmsg -ptr byte' -float0 -ew
endif


# ----------------------------------------------------------------------------
# Create dependency list.
# ----------------------------------------------------------------------------

INITIAL    = initial

OPERATORS  = polyops	operators

FFT	   = temfftd	temffts		fftpack		fourier

MEMORY     = mapping	family          message

MESH       = RCM

ALL        = $(INITIAL) $(OPERATORS) $(FFT) $(MEMORY) $(MESH)
     
LMD        = $(foreach routine, $(ALL), $(LIB)($(routine).o))

# ----------------------------------------------------------------------------

# -- Make library, default action.

$(LIB): $(LMD)
ifeq ($(RANLIB),ranlib)
	$(RANLIB) $(LIB)
endif

# -- Install library & associated header files.

install:
	if test ! -d $(SEMHOME)/lib 		;		then	\
		mkdir $(SEMHOME)/lib 		; 		fi;	\
	if test ! -d $(SEMHOME)/lib/$(ARCH) 	;		then	\
		mkdir $(SEMHOME)/lib/$(ARCH)	; 		fi;	\
	if test -r $(LIB) 			;		then	\
		rm -f $(SEMHOME)/lib/$(ARCH)/$(LIB)	;		\
		cp $(LIB) $(SEMHOME)/lib/$(ARCH)	;	fi

	rm -f $(SEMHOME)/include/femdef.h;				\
	rm -f $(SEMHOME)/include/femlib.h;				\
	rm -f $(SEMHOME)/include/Femlib.h;				\

	cp femdef.h    $(SEMHOME)/include	;			\
	cp femlib.h    $(SEMHOME)/include	;			\
	cp Femlib.h    $(SEMHOME)/include	;			\

initial.o: defaults.h

testpoly: testpoly.c
	cc $(CFLAGS) -o testpoly testpoly.c -I$(SEMHOME)/include -I. \
	-L$(SEMHOME)/lib -lfem -lalp -lm -lmalloc

testFFT: testFFT.C
	$(CXX) -o $@ testFFT.C temfft.o -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lftn -lm

tar:
	tar cvf femlib.tar *

clean:
	rm -rf *.o *.a *~ ii_files *.F

# -- Special compiles for .F FFT routines, there is a problem otherwise.

ifeq ($(ARCH),IRIX)
temfftd.o: temfftd.F
	f77 -c -g temfftd.F
temffts.o: temffts.F
	f77 -c -O2 temffts.F
endif

ifeq ($(ARCH),IRIX64)
temfftd.o: temfftd.F
	f77 -c -O2 temfftd.F
temffts.o: temffts.F
	f77 -c -O2 temffts.F
endif

ifeq ($(ARCH),OSF1)
temfftd.o: temfftd.F
	f90 -c -fast -inline speed -tune host temfftd.F -O5 -unroll 2
temffts.o: temffts.F
	f90 -c -fast -inline speed -tune host temffts.F -O5 -unroll 2
endif

ifeq ($(ARCH),Linux)
temfftd.o: temfftd.F
	g77 -c -O3 temfftd.F
temffts.o: temffts.F
	g77 -c -O3 temffts.F
endif

ifeq ($(ARCH),UNIX_System_V)
temfftd.o: temfftd.F
	frt -c -Wv,-m3,-te temfftd.F
temffts.o: temffts.F
	frt -c -Wv,-m3,-te temffts.F
endif

ifeq ($(ARCH),SUPER-UX)
temfftd.o: temfftd.F
	f77 -c -Wf'-pvctl,fullmsg -ptr byte' -float0 -ew temfftd.F
temffts.o: temffts.F
	f77 -c -Wf'-pvctl,fullmsg -ptr byte' -float0 -ew temffts.F
endif
