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

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)			# -- SGI 32-bit systems.
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -g
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -L$(SEM)/lib/$(ARCH) -lfem -lalp \
	   $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)			# -- SGI 64-bit systems.
CXX	 = CC
LD	 = CC
ifdef DEBUG
  OPT    = -g
else
  OPT    =  -r10000 -Ofast=IP30 -OPT:roundoff=3
endif
CXXFLAGS = $(DEFINES) -mips4 -64 $(OPT)
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -64 -mips4 $(OPT) -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),OSF1)			# -- DEC OSF1 Alpha.
CXX	 = cxx
LD	 = cxx
ifdef DEBUG
  OPT    = -g
else
  OPT    = -tune host -O2 -Olimit=1000
endif
CXXFLAGS = $(DEFINES) $(OPT)
F77LIBS  = -ldxml -lfor
ifdef MPI
  CXXFLAGS += -thread_safe
  LDFLAGS   = $(CXXFLAGS) -threads -L/usr/local/mpi/lib/alpha/ch_shmem	\
	     -L$(SEM)/lib/$(ARCH) -lfem_mp -lalp_mp $(F77LIBS) -lm -lmpi
else
  LDFLAGS  = $(CXXFLAGS) -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif
endif

ifeq ($(ARCH),Linux-alpha)		# -- Linux alpha.
CXX      = g++
LD       = g++
F77LIBS  = -lcxml -lots -lfor -lg2c
#F77LIBS = -llapack -lgemm -lblas -lots -lfor -lg2c
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3
endif
CXXFLAGS = $(DEFINES) $(OPT)
ifdef MPI
  LDFLAGS = -L$(SEM)/lib/$(ARCH) -L/usr/lib -L$(MPIHOME)/lib  \
	    -lfem_mp -lalp $(F77LIBS) -lcpml 		      \
            -lmpi -ltrillium -ltstdio -lt -largs
#           -lmpich
else
  LDFLAGS = -L$(SEM)/lib/$(ARCH) -L/usr/local/lib    \
            -lfem -lalp $(F77LIBS) -lcpml
endif
endif

ifeq ($(ARCH),Linux-i586)		# -- Linux Intel.
CXX	 = g++
LD	 = g++
F77LIBS  = -llapack -lblas -lg2c
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3
endif
CXXFLAGS = $(DEFINES) $(OPT)
LDFLAGS  = $(OPT) -L$(SEM)/lib/$(ARCH) -L/usr/local/lib	\
            -lfem -lalp $(F77LIBS)
endif


ifeq ($(ARCH),UNIX_System_V)    	# -- Fujitsu VPP300.
CC       = vcc
CXX      = /home/home01/565/hmb565/C++/bin/CC
LD       = /home/home01/565/hmb565/C++/bin/CC
OPTIM    = -Wv,-m3,-Om,-te -Kpopt
CXXFLAGS = $(DEFINES)  $(OPTIM)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS)
ifdef MPI
  LIBDIR   = -L$(SEM)/lib/$(ARCH)
  LIBDIR  += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
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

ifeq ($(ARCH),SUPER-UX)         	# -- NEC SX-4.
DEFINES += -D_BUILTIN_
LIBDIR   = -L$(SEM)/lib/$(ARCH) 
CC       = cc
CXX      = CC
LD       = CC
CFLAGS   = $(DEFINES) -pvctl,fullmsg -hacct -Onooverlap        \
           -hfloat0 -hnostkchk
CXXFLAGS = $(CFLAGS) $(OPTIM) +a1 -ispeed
LDFLAGS  = $(CFLAGS) $(CPPFLAGS) $(LIBDIR) $(LIBS)
ifdef MPI
  LIBS   = -lfem_mp -llapack_64 -lalp_mp -lblas_64 -lfft_64	\
           -lu90sxe -lp90sx -lv90sxe -lf90sxe -li90sx -lm -lmpi -lpthread
else
  LIBS   = -lfem -llapack_64 -lalp -lblas_64 -lfft_64	\
           -lu90sxe -lp90sx -lv90sxe -lf90sxe -li90sx -lm
endif
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
NS_OBJ = drive.o NS.o dnsanalysis.o

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
