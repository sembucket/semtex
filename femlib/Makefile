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
CFLAGS     = $(DEFINES) -mips4 -64 -r10000 -Ofast=IP30
FFLAGS     = $(DEFINES) -mips4 -64 -r10000 -Ofast=IP30 -OPT:roundoff=3
endif

ifeq ($(ARCH),OSF1)			# -- DEC Alpha OSF1
FC	   = f90
CXX	   = cxx
CC         = cc
CFLAGS     = $(DEFINES) -fast -inline speed -tune host -O4
FFLAGS     = $(DEFINES) -fast -inline speed -tune host -O5 -unroll 2
ifdef MPI
  LIB       = libfem_mp.a
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
CFLAGS     = $(DEFINES) -Wv,-m3,-Om,-te
FFLAGS     = $(DEFINES) -Wv,-m3,-te
ifdef MPI
  LIB       = libfem_mp.a
  DEFINES  += -DMPI
  CPPFLAGS += -I/usr/lang/mpi/include
endif
endif

ifeq ($(ARCH),SUPER-UX)                 # -- NEC SX-4
DEFINES   += -D_BUILTIN_
CFLAGS     = $(DEFINES) -pvctl,fullmsg -hfloat0 -hnostkchk
FFLAGS     = $(DEFINES) -Cvopt -Wf'-pvctl,fullmsg -ptr byte' -float0 -ew
ifdef MPI
  LIB       = libfem_mp.a
  DEFINES  += -DMPI
  CFLAGS   += 
  FFLAGS   += -P multi
endif
endif


# ----------------------------------------------------------------------------
# Create dependency list.
# ----------------------------------------------------------------------------

INITIAL    = initial

OPERATORS  = polyops	operators	matops

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
	rm -rf *.O *.o *.a *~ ii_files *.F

# -- Special machine-dependent optimizations.

ifeq ($(ARCH),OSF1)
temfftd.o: temfftd.F
	f90 -c -fast -inline speed -tune host temfftd.F -O5 -unroll 2
temffts.o: temffts.F
	f90 -c -fast -inline speed -tune host temffts.F -O5 -unroll 2
endif

# -- Rules to help with .F routines.

matops.o:  matops.F
temfftd.o: temfftd.F
temffts.o: temffts.F
