##############################################################################
# Makefile for dns, spectral element solver for unsteady Navier--Stokes
# problems.
#
# $Id$
##############################################################################

SEM   = ..
VPATH = $(SEM)/src
ARCH  = $(shell uname -s)
MACH  = $(shell uname -m)
ifeq ($(ARCH),Linux)
  ifeq ($(MACH),alpha)
    ARCH = Linux-alpha
  endif
  ifeq ($(MACH),i586)
    ARCH = Linux-i586
  endif
endif

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
endif
ifdef STOKES
  DEFINES += -DSTOKES
endif
ifdef MPI
  MPIHOME = /usr/local/lam
  PROG    = dns_mp
else
  PROG = dns
endif

CPPFLAGS = -I. -I$(SEM)/include
LIBDIR   = -L$(SEM)/lib/$(ARCH)

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),Linux-alpha)              # -- Linux alpha.
DEFINES += -D_Linux_alpha
CXX      = g++
LD       = g++
LIBDIR  += /usr/local/lib
LDFLAGS  = $(LIBDIR)
ifdef STATIC
  F77LIBS  = -lcxml -lfor -lFutil -lUfor -lots -lg2c
  LDFLAGS + -static
else
  F77LIBS  = -lcxml -lfor -lots -lg2c
endif
ifdef NOCXML     # -- For testing, avoid Compaq math libs.
  F77LIBS =  -llapack -lgemm -lblas -lots -lfor -lg2c  
endif
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3
endif
CXXFLAGS = $(DEFINES) $(OPT)
ifdef MPI
  LIBDIR  += -L$(MPIHOME)/lib
  LDFLAGS += -lfem_mp -lalp $(F77LIBS) -lcpml \
             -lmpi -ltrillium -ltstdio -lt -largs
#           -lmpich
else
  LDFLAGS  += -lfem -lalp $(F77LIBS) -lcpml
endif
endif

ifeq ($(ARCH),Linux-i586)               # -- Linux Intel.
CXX      = g++
LD       = g++
F77LIBS  = -llapack -lblas -lg2c
LIBDIR  += -L/usr/local/lib
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3
endif
CXXFLAGS = $(DEFINES) $(OPT)
LDFLAGS  = $(OPT) $(LIBDIR) -lfem -lalp $(F77LIBS)
endif

ifeq ($(ARCH),IRIX)		# -- SGI 32-bit systems.
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -fullwarn -g
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = $(LIBDIR) -lfem -lalp $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI 64-bit systems.
CXX	 = CC
LD	 = CC
ifdef DEBUG
  OPT    = -g
else
  OPT    =  -r10000 -Ofast=IP30 -OPT:roundoff=3
endif
CXXFLAGS = $(DEFINES) -mips4 -64 $(OPT)
F77LIBS  = -lcomplib.sgimath -lftn
LDFLAGS  = -64 -mips4 $(OPT) $(LIBDIR) -lfem -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha.
CXX	 = cxx
LD	 = cxx
CXXFLAGS =  $(DEFINES) -tune host -O2
F77LIBS  = -ldxml -lfor
ifdef MPI
  CXXFLAGS += -thread_safe
  LIBDIR   += -L/usr/local/mpi/lib/alpha/ch_shmem
  LDFLAGS   = $(CXXFLAGS) -threads  $(LIBDIR) \
	      -lfem_mp -lalp_mp $(F77LIBS) -lm -lmpi
else
  LDFLAGS  = $(CXXFLAGS) $(LIBDIR) -lfem -lalp $(F77LIBS) -lm
endif
endif

ifeq ($(ARCH),UNIX_System_V)    # -- Fujitsu VPP300.
CC       = vcc
CXX      = /usr/uxp/C++/bin/CC
LD       = /usr/uxp/C++/bin/CC
OPTIM    = -Wv,-m3,-Om,-te -Kpopt
CXXFLAGS = $(DEFINES)  $(OPTIM)
ifdef MPI
  LIBDIR  += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
  LIBDIR  += -L/usr/lang/mpi/lib
  LIBS     = -lfem_mp -lalp
  LIBS    += -llapack -lblasvpp  -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
  LIBS    += -lmpi -lmp
  OPTIM   += -Wl,-P
else
  LIBDIR  += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
  LIBS     = -lfem -lalp
  LIBS    += -llapack -lblasvpp  -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
endif
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS)
endif

ifeq ($(ARCH),SUPER-UX)         # -- NEC SX-4.
DEFINES += -D_BUILTIN_
CC       = cc
CXX      = CC
LD       = CC
DEFINES += -D_SX
CFLAGS   = $(DEFINES) -pvctl,fullmsg -hacct -Onooverlap -hfloat0 -hnostkchk
CXXFLAGS = $(CFLAGS) $(OPTIM) +a1 -ispeed
ifdef MPI
  LIBS   = -lfem_mp -lalp_mp -llapack_64 -lblas_64 -lfft_64	\
           -lu77sxe -lp77sx -lv77sxe -lf77sxe -li77sx -lm -lmpi -lpthread
else
  LIBS   = -lfem -lalp -llapack_64 -lblas_64 -lfft_64	\
           -lu77sxe -lp77sx -lv77sxe -lf77sxe -li77sx -lm
endif
LDFLAGS  = $(CFLAGS) $(CPPFLAGS) $(LIBDIR) $(LIBS)
endif

# ----------------------------------------------------------------------------
#
SEMFILES = analysis auxfield BCmgr boundary bsys condition domain   \
	   element feml field geometry history integration  	\
	   matrix mesh misc particle pressure statistics
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = Sem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
NS_FILES = drive NS dnsanalysis
NS_OBJ   = $(addsuffix .o,$(NS_FILES))

$(PROG): sem $(SEMOBJ) $(NS_OBJ) 
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
	$(CXX) -I$(SEM)/include -c -hnostkchk +a1 -hacct -hnovector -p \
	$(VPATH)/mesh.C
BCmgr.o: $(VPATH)/BCmgr.C
	$(CXX) -I$(SEM)/include -c -hnostkchk +a1 -hacct -hnovector -p \
	$(VPATH)/BCmgr.C
bsys.o: $(VPATH)/bsys.C
	$(CXX) -I$(SEM)/include -c -hnostkchk +a1 -hacct -hnovector -p \
	$(VPATH)/bsys.C
condition.o: $(VPATH)/condition.C
	$(CXX) -I$(SEM)/include -c -hnostkchk +a1 -hacct -hnovector -p \
	$(VPATH)/condition.C
domain.o: $(VPATH)/domain.C
	$(CXX) -I$(SEM)/include -c -hnostkchk +a1 -hacct -hnovector -p \
	$(VPATH)/domain.C
dnsanalysis.o: dnsanalysis.C
	$(CXX) -I$(SEM)/include -I. -c -hnostkchk +a1 -hacct -hnovector -p \
	dnsanalysis.C
endif

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
	rm -rf ptrepository
