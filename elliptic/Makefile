##############################################################################
# Makefile for spectral element Helmholtz solver.
#
# $Id$
##############################################################################

SEM      = ..
VPATH    = $(SEM)/src

ARCH     = $(shell uname -s)

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
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
CXXFLAGS =  $(DEFINES) -64 -mips4 -g
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -g -64 -mips4 -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CXX	 = cxx
LD	 = cxx
CXXFLAGS =  $(DEFINES) -migrate -tune host -g -w0
F77LIBS  = -ldxml
ifdef MPI
  CXXFLAGS += -thread_safe
  LDFLAGS  = $(CXXFLAGS) -threads -L/home/hmb/mpi/lib/alpha/ch_shmem	\
             -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm -lmpi
else
  LDFLAGS  = $(CXXFLAGS) -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CXX	 = g++
LD	 = g++
CXXFLAGS = $(DEFINES) -gstabs+
F77LIBS  = -llapack -lblas -lf2c
LDFLAGS  = -gstabs+ -L$(SEM)/lib/$(ARCH) -L/usr/local/lib	\
	   -lfem -lalp $(F77LIBS)
endif

ifeq ($(ARCH),UNIX_System_V)    # -- Fujitsu VPP300.
LIBDIR   = -L$(SEM)/lib/$(ARCH) 
LIBDIR  += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
LIBS     = -lfem -lalp
LIBS    += -llapack -lblasvpp  -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
CC       = vcc
CXX      = /usr/uxp/C++/bin/CC
LD       = /usr/uxp/C++/bin/CC
OPTIM    = -Wv,-m3,-Ps -Kpopt
#OPTIM    = -Wv,-ad
CXXFLAGS = $(DEFINES)  $(OPTIM)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS)
endif

ifeq ($(ARCH),SUPER-UX)         # -- NEC SX-4.
DEFINES += -D_BUILTIN_
LIBDIR   = -L$(SEM)/lib/$(ARCH) 
CC       = cc
CXX      = CC
LD       = CC
CFLAGS   = $(DEFINES) -pvctl,fullmsg -hacct -Onooverlap        \
           -hfloat0 -hnostkchk
CXXFLAGS = $(CFLAGS) $(OPTIM) +a1
LDFLAGS  = $(CFLAGS) $(CPPFLAGS) $(LIBDIR) $(LIBS)
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

SEMFILES  = auxfield BCmgr boundary domain element feml	\
	    field geometry integration matrix mesh misc pressure
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h


# ----------------------------------------------------------------------------
#
ELLOBJ = drive.o helmholtz.o

elliptic: sem $(SEMHDR) $(SEMOBJ) $(ELLOBJ)
	$(LD) -o elliptic $(ELLOBJ) $(SEMOBJ) $(CPPFLAGS) $(LDFLAGS)

# ----------------------------------------------------------------------------
#
$(SEMOBJ): $(SEMHDR)

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
condition.o: $(VPATH)/condition.C
        $(CXX) -I$(SEM)/include -c -hnostkchk +a1 -hacct -hnovector -p \
        $(VPATH)/condition.C
endif

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
	rm -rf ptrepository


