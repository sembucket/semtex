##############################################################################
# Makefile for ns, spectral element solver for unsteady Navier--Stokes
# problems.
#
# $Id$
##############################################################################

SEM   = ..
VPATH = $(SEM)/src
ARCH  = $(shell uname -s)

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
endif
ifdef STOKES
  DEFINES += -DSTOKES
endif

CPPFLAGS = -I. -I$(SEM)/include

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI 32-bit systems.
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -g
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -L$(SEM)/lib/$(ARCH) -lfem -lalp \
	   $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI 64-bit systems.
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -mips4 -64 -Ofast=IP30 -r10000
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -64 -mips4 -Ofast=IP30 -r10000 	\
           -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha.
CXX	 = cxx
LD	 = cxx
CXXFLAGS =  $(DEFINES) -tune host -O2
F77LIBS  = /usr/lib/libdxml.a -lfor
ifdef MPI
  CXXFLAGS += -thread_safe
  LDFLAGS   = $(CXXFLAGS) -threads -L/usr/local/mpi/lib/alpha/ch_shmem	\
	     -L$(SEM)/lib/$(ARCH) -lfem_mp -lalp_mp $(F77LIBS) -lm -lmpi
else
  LDFLAGS  = $(CXXFLAGS) -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CXX	 = g++
LD	 = g++
CXXFLAGS = $(DEFINES) -O3
F77LIBS  = -llapack -lblas -lf2c
LDFLAGS  = -gstabs+ -L$(SEM)/lib/$(ARCH) -L/usr/local/lib	\
	   -lfem -lalp $(F77LIBS)
endif

ifeq ($(ARCH),UNIX_System_V)    # -- Fujitsu VPP300.
CC       = vcc
CXX      = /usr/uxp/C++/bin/CC
LD       = /usr/uxp/C++/bin/CC
OPTIM    = -Wv,-m3,-Om,-te
CXXFLAGS = $(DEFINES)  $(OPTIM)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS)
ifdef MPI
  LIBDIR   = -L$(SEM)/lib/$(ARCH)
  LIBDIR  += -L/usr/lang/mpi/lib
  LIBS     = -lfem_mp -lalp
  LIBS    += -llapack -lblasvpp  -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
  LIBS    += -lmpi -lmp
  LDFLAGS += -Wl,-P
else
  LIBDIR   = -L$(SEM)/lib/$(ARCH) 
  LIBDIR  += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
  LIBS     = -lfem -lalp
  LIBS    += -llapack -lblasvpp  -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
endif
endif

ifeq ($(ARCH),SUPER-UX)         # -- NEC SX-4.
LIBDIR   = -L$(SEM)/lib/$(ARCH) 
CC       = cc
CXX      = CC
LD       = CC
OPTIM    = -pvctl,fullmsg -hacct -hvector -p -ptn -Onooverlap
CXXFLAGS = $(DEFINES)  $(OPTIM)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS) -hacct
ifdef MPI
  LIBS   = -lfem_mp -llapack_64 -lalp_mp -lblas_64 -lfft_64	\
           -lu77sxe -lp77sx -lv77sxe -lf77sxe -li77sx -lm -lmpi -lpthread
else
  LIBS   = -lfem -llapack_64 -lalp -lblas_64 -lfft_64	\
           -lu77sxe -lp77sx -lv77sxe -lf77sxe -li77sx -lm
endif
endif

# ----------------------------------------------------------------------------
#
SEMFILES  = analysis auxfield BCmgr boundary condition domain 	\
	    element feml field geometry history integration  	\
	    matrix mesh misc particle pressure statistics
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
NS_OBJ = drive.o NS.o

dns:	sem $(SEMOBJ) $(NS_OBJ) 
	$(LD) -o $@ $(NS_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
#
$(SEMOBJ): $(SEMHDR)

#-----------------------------------------------------------------------------
sem:
	cd $(VPATH); gmake install

# ----------------------------------------------------------------------------
# Special machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	CC -I$(SEM)/include -c -pvctl,fullmsg -hacct -hvector -p -ptn   \
	$(VPATH)/mesh.C
BCmgr.o: $(VPATH)/BCmgr.C
	CC -I$(SEM)/include -c -pvctl,fullmsg -hacct -hvector -p -ptn   \
	$(VPATH)/BCmgr.C
endif

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
	rm -rf ptrepository
