##############################################################################
# Makefile for spectral element Helmholtz solver.
#
# $Id$
##############################################################################

SEM   = ..
VPATH = $(SEM)/src

ARCH = $(shell uname -s)
MACH = $(shell uname -m)
ifeq ($(ARCH),Linux)
  ifeq ($(MACH),alpha)
    ARCH = Linux-alpha
  endif
  ifeq ($(MACH),i586)
    ARCH = Linux-i586
  endif
endif

ifdef DEBUG
  DEFINES = -DDEBUG
else
  DEFINES =
endif
ifdef MPI
  MPIHOME = /usr/local/lam
  PROG    = elliptic_mp
else
  PROG = elliptic
endif

CPPFLAGS = -I. -I$(SEM)/include
LIBDIR   = -L$(SEM)/lib/$(ARCH)

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI 32-bit systems.
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -g
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = $(LIBDIR) -lfem -lalp $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI 64-bit systems.
CXX	 = CC
LD	 = CC
CXXFLAGS =  $(DEFINES) -64 -g
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -g -64 $(LIBDIR) -lfem -lalp $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CXX	 = cxx
LD	 = cxx
CXXFLAGS = $(DEFINES) -migrate -tune host -Olimit 1000 -g -w0
F77LIBS  = -ldxml
ifdef MPI
  CXXFLAGS += -thread_safe
  LDFLAGS  = $(CXXFLAGS) -threads -L/usr/local/mpi/lib/alpha/ch_shmem	\
             $(LIBDIR) -lfem_mp -lalp_mp $(F77LIBS) -lm -lmpi
else
  LDFLAGS  = $(CXXFLAGS) $(LIBDIR) -lfem -lalp $(F77LIBS) -lm
endif
endif

ifeq ($(ARCH),Linux-i586)		# -- Linux Intel system.
CXX = g++
LD  = g++
F77LIBS = -llapack -lblas -lg2c
ifdef DEBUG
  OPT = -gstabs+
else
  OPT = -O3
endif	
CXXFLAGS = $(DEFINES) $(OPT)
LDFLAGS  = $(OPT) $(LIBDIR) -L/usr/local/lib -lfem -lalp $(F77LIBS)
endif

ifeq ($(ARCH),Linux-alpha)		# -- Linux Alpha system.
CXX     = g++
LD      = g++
F77LIBS = -lcxml -lots -lfor -lg2c
LIBDIR += -L/usr/local/lib
ifdef DEBUG
  OPT = -gstabs+
else
  OPT = -O3
endif
CXXFLAGS = $(DEFINES) $(OPT)
ifdef MPI
  LDFLAGS = $(LIBDIR) -L/usr/lib -L$(MPIHOME)/lib \
            -lfem_mp -lalp $(F77LIBS) -lcpml                 \
	    -lmpi -ltrillium -ltstdio -lt -largs 
#           -lmpich
else
  LDFLAGS = $(LIBDIR) -lfem -lalp $(F77LIBS) -lcpml
endif
endif

ifeq ($(ARCH),UNIX_System_V)    # -- Fujitsu VPP300.
LIBDIR  += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
LIBS     = -lfem -lalp
LIBS    += -llapack -lblasvpp  -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
CC       = vcc
CXX      = /home/home01/565/hmb565/C++/bin/CC
LD       = $(CXX)
OPTIM    = -Wv,-m3,-Ps -Kpopt
#OPTIM    = -Wv,-ad
CXXFLAGS = $(DEFINES)  $(OPTIM)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS)
endif

ifeq ($(ARCH),SUPER-UX)         # -- NEC SX-4.
CXX       = c++
LD        = c++ -verbose -V -p
DEFINES  += -D_BUILTIN_ -D_SX
CPPFLAGS += -I/usr/include/CC
OPT       = -O nooverlap -Nstkchk
NOVEC     = -c -Nstkchk -acct -dir novec
CXXFLAGS  = -V $(DEFINES) -pvctl,loopcnt=20000,fullmsg -acct $(OPT)
FORTLIBS  = -llapack_64 -lblas_64 -lfft_64
ifdef MPI
  LIBS    = -lfem_mp -lalp_mp $(FORTLIBS) -lm -lC -lmpi -lpthread
else
  LIBS    = -lfem -lalp $(FORTLIBS) -lm -lC
endif
LDFLAGS   = -verbose -V -p -f90libew $(CXXFLAGS) $(CPPFLAGS) $(LIBDIR) $(LIBS)
endif

# ----------------------------------------------------------------------------
#

SEMFILES  = auxfield BCmgr boundary bsys condition domain element feml	\
	    field geometry integration matrix mesh misc pressure  
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h


# ----------------------------------------------------------------------------
#
ELLOBJ = drive.o helmholtz.o

elliptic: sem $(SEMHDR) $(SEMOBJ) $(ELLOBJ)
	$(LD) -o $(PROG) $(ELLOBJ) $(SEMOBJ) $(CPPFLAGS) $(LDFLAGS)

# ----------------------------------------------------------------------------
#
$(SEMOBJ): $(SEMHDR)

sem:
	cd $(VPATH); gmake install

# ----------------------------------------------------------------------------
# Special machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
	rm -rf ptrepository
