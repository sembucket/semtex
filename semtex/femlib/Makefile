##############################################################################
# Remake finite element routine library (libfem.a).  Gnu make required.
#
# $Id$
##############################################################################

ARCH = $(shell uname -s)

OPTIONS  = -g

FEMHOME  = ..
LIB      = libfem.a

DEFINES  =
ifdef DEBUG
  DEFINES += -DDEBUG
endif

CC	 = cc
CXX	 = CC
FC	 = f77
CFLAGS	 = $(OPTIONS)
FFLAGS   = $(OPTIONS)
ARFLAGS	 = cr
CPPFLAGS = -I$(FEMHOME)/include
RANLIB   =

# -------------- Compiler options for specific machines. ---------------------

ifeq ($(ARCH),IRIX)			# -- SGI IRIX 32-bit
CFLAGS     = $(DEFINES) -fullwarn -O2
FFLAGS     = $(DEFINES) -O2
F77LIBS    = -lcomplib.sgimath -lftn
endif

ifeq ($(ARCH),IRIX64)			# -- SGI IRIX 64-bit
CFLAGS     = $(DEFINES) -fullwarn -mips4 -O2
FFLAGS     = $(DEFINES) -mips4 -O2
F77LIBS    = -lcomplib.sgimath -lftn
endif

ifeq ($(ARCH),OSF1)			# -- DEC Alpha OSF1
FC	   = f90
CXX	   = cxx
CFLAGS     = $(DEFINES) -migrate -fast -tune host -O4
FFLAGS     = $(DEFINES) -migrate -fast -tune host -O5
F77LIBS    = -ldxml
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
LIBDIR     = -L$(SEM)/lib/$(ARCH) 
LIBDIR    += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib 
LIBS       = -lfem -lalp
LIBS      += -lblasvpp -lm -lfj90fv -lfj90 -lfj90f -ljsp -lvfl
CC         = vcc
FC         = frt
CXX        = /usr/uxp/C++/bin/CC
LD         = /usr/uxp/C++/bin/CC
OPTIM      = -Wv,-m3,-Om,-te
CXXFLAGS   = $(OPTIM)
CFLAGS     = $(OPTIM)
LDFLAGS    = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS)
FFLAGS     = -Wv,-m3,-te
RANLIB     = 
endif

ifeq ($(ARCH),SUPER-UX)                 # -- NEC SX-4
CC         = cc
FC         = f77
CFLAGS     = -pvctl,fullmsg -O nooverlap -float0
FFLAGS     = -Wf'-pvctl,fullmsg -ptr byte' -float0 -ew
RANLIB     = 
endif


# ----------------------------------------------------------------------------
# Create dependency list.
# ----------------------------------------------------------------------------

INITIAL    = initial

OPERATORS  = polyops	operators

FFT	   = temfftd	temffts		fftpack		fourier

MEMORY     = mapping	family

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
	if test ! -d $(FEMHOME)/lib 		;		then	\
		mkdir $(FEMHOME)/lib 		; 		fi;	\
	if test ! -d $(FEMHOME)/lib/$(ARCH) 	;		then	\
		mkdir $(FEMHOME)/lib/$(ARCH)	; 		fi;	\
	if test -r $(LIB) 			;		then	\
		rm -f $(FEMHOME)/lib/$(ARCH)/$(LIB)	;		\
		cp $(LIB) $(FEMHOME)/lib/$(ARCH)	;	fi

	rm -f $(FEMHOME)/include/femdef.h;				\
	rm -f $(FEMHOME)/include/femlib.h;				\
	rm -f $(FEMHOME)/include/Femlib.h;				\

	cp femdef.h    $(FEMHOME)/include	;			\
	cp femlib.h    $(FEMHOME)/include	;			\
	cp Femlib.h    $(FEMHOME)/include	;			\

initial.o: defaults.h

testpoly: testpoly.c
	cc $(CFLAGS) -o testpoly testpoly.c -I$(FEMHOME)/include -I. \
	-L$(FEMHOME)/lib -lfem -lalp -lm -lmalloc

testFFT: testFFT.C
	$(CXX) -o $@ testFFT.C temfft.o -I$(FEMHOME)/include -I. -g \
	-L$(FEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lftn -lm

tar:
	tar cvf femlib.tar *

clean:
	rm -rf *.o *.a *~ ii_files *.F

# -- Special compiles for .F FFT routines, there is a problem otherwise.

ifeq ($(ARCH),IRIX)
temfftd.o: temfftd.F
	f77 -c -O2 temfftd.F
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
	f90 -c -migrate -fast -tune host temfftd.F -O5
temffts.o: temffts.F
	f90 -c -migrate -fast -tune host temffts.F -O5
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
        f77 -c -Wf'-pvctl,fullmsg' -float0 -ew temfftd.F
temffts.o: temffts.F
        f77 -c -Wf'-pvctl,fullmsg' -float0 -ew temffts.F
endif
