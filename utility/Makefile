##############################################################################
# Makefile for spectral element solver utilities.
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
  DEFINES  = -DDEBUG
else
  DEFINES  =
endif

CPPFLAGS = -I. -I$(SEM)/include
LIBDIR   = -L$(SEM)/lib/$(ARCH) 

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x
CC       = cc
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -fullwarn -mips2 -g
FORTLIBS = -lcomplib.sgimath -lftn
LDFLAGS  = -lfpe -mips2 $(LIBDIR) -lfem -lalp $(FORTLIBS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CC       = cc
CXX	 = CC
LD	 = CC
ifdef DEBUG
  OPT    = -mips4 -64 -g
else
  OPT    = -mips4 -64 -Ofast=IP30 -r10000
endif
CXXFLAGS = $(DEFINES) $(OPT)
FORTLIBS = -lcomplib.sgimath -lftn
LDFLAGS  = $(OPT) $(LIBDIR) -lfem -lalp $(FORTLIBS) -lmalloc -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CC       = cc
CXX	 = cxx
LD	 = cxx
CXXFLAGS =  $(DEFINES) -migrate -tune host -g
FORTLIBS = -ldxml
LDFLAGS  = $(CXXFLAGS) $(LIBDIR) -lfem -lalp $(FORTLIBS) -lm
endif

ifeq ($(ARCH),Linux-i586)	# -- Linux Intel system.
CC	 = gcc
CXX	 = g++
LD	 = g++
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O2
endif
CXXFLAGS = $(DEFINES) $(OPT)
CFLAGS   = $(DEFINES) $(OPT)
FORTLIBS = -llapack -lblas -lg2c
LDFLAGS  = $(LIBDIR) -L/usr/local/lib -lfem -lalp $(FORTLIBS) -lm
endif

ifeq ($(ARCH),Linux-alpha)	# -- Linux alpha system.
CC	 = ccc
CXX	 = g++
LD	 = g++
ifdef DEBUG
  OPT    = -g
else
  OPT    = -O2
endif
CXXFLAGS = $(DEFINES) $(OPT)
FORTLIBS = -lcxml -lots -lfor -lg2c
LDFLAGS  = $(LIBDIR) -L/usr/local/lib -lfem -lalp $(FORTLIBS) -lm
endif

ifeq ($(ARCH),UNIX_System_V)    # -- Fujitsu VPP300.
LIBDIR   = -L$(SEM)/lib/$(ARCH) 
LIBDIR  += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
LIBS     = -lfem -lalp
LIBS    += -llapack -lblasvpp -lm -lfj90fv -lfj90 -lfj90f -ljsp -lvfl
CC       = cc
CXX      = /home/home01/565/hmb565/C++/bin/CC
LD       = $(CXX)
ifdef DEBUG
  OPT    = -g
else
  OPT    = -Wv,-m3,-Ps 
endif
CXXFLAGS = $(DEFINES)  $(OPT)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS)
endif

ifeq ($(ARCH),SUPER-UX)         # -- NEC SX-4.
CXX       = c++
LD        = c++
DEFINES  += -D_BUILTIN_ -D_SX
CPPFLAGS += -I/usr/include/CC
OPT       = -O nooverlap -Nstkchk
NOVEC     = -c -Nstkchk -acct -dir novec
CXXFLAGS  = -V $(DEFINES) -pvctl,loopcnt=20000,fullmsg -acct $(OPT)
FORTLIBS  = -llapack_64 -lblas_64 -lfft_64
ifdef MPI
  LIBS    = -lfem_mp -lalp $(FORTLIBS) -lm -lC -lmpi -lpthread
else
  LIBS    = -lfem -lalp $(FORTLIBS) -lm -lC
endif
LDFLAGS   =  -verbose -V -p -f90libew $(CXXFLAGS) $(CPPFLAGS) $(LIBDIR) $(LIBS)
endif

# ----------------------------------------------------------------------------
headers:
	cd $(VPATH); make install

all: meshpr  enumerate calc   addfield  sem2nek  sem2tec \
     convert compare   noiz   project   interp   probe   \
     rstress rectmesh  transform

# ----------------------------------------------------------------------------
# Build meshpr utility.
MPROBJ = feml.o mesh.o meshpr.o

meshpr:	headers $(MPROBJ)
	$(LD) -o $@ $(MPROBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build enumerate utility.
ENUOBJ = feml.o mesh.o enumerate.o

enumerate:	headers $(ENUOBJ)
	$(LD) -o $@ $(ENUOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build calc utility.

calc: calc.o
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build addvort utility.
ADDOBJ = addvort.o auxfield.o BCmgr.o boundary.o bsys.o condition.o 	\
	 domain.o element.o feml.o field.o geometry.o integration.o	\
	 mesh.o misc.o pressure.o

addvort:	$(ADDOBJ)
	$(LD) -o $@ $(ADDOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build addvortr utility.
ADROBJ = addvortr.o auxfield.o BCmgr.o boundary.o domain.o	\
	 element.o feml.o field.o geometry.o  integration.o mesh.o	\
	 misc.o	pressure.o \

addvortr:	$(ADROBJ)
	$(LD) -o $@ $(ADROBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build addfield utility.
ADFOBJ = addfield.o auxfield.o BCmgr.o boundary.o bsys.o condition.o 	\
	 domain.o element.o feml.o field.o geometry.o integration.o	\
	 mesh.o misc.o pressure.o

addfield:	$(ADFOBJ)
	$(LD) -o $@ $(ADFOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build semtex --> NEKTON input file converter.
S2NOBJ = feml.o mesh.o sem2nek.o 

sem2nek: $(S2NOBJ)
	$(LD) -o $@ $(S2NOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build convert, binary <--> ASCII field dump converter.

convert: convert.c
	$(CC) -o $@ $@.c

# ----------------------------------------------------------------------------
# Build sem2tec, field dump --> tecplot input file converter.

sem2tec: sem2tec.o
	$(CC) $(CPPFLAGS) -o $@ $@.c $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build compare, exact solution comparison/restart utility.
COMOBJ = compare.o bsys.o BCmgr.o boundary.o condition.o feml.o mesh.o	\
	 element.o auxfield.o geometry.o pressure.o integration.o misc.o

compare: $(COMOBJ)
	$(LD) -o $@ $(COMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build noiz, velocity field perturber.

noiz: noiz.c
	$(CC) $(CPPFLAGS) -o $@ $@.c $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build moden, to compute TKE distribution for a given mode.

moden: moden.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.c $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build project, field dump projection to higher/lower orders.

project: project.C
	$(CXX) $(CPPFLAGS) -g -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build transform, which does forward/inverse Fourier/Legendre transforms.

transform: transform.C
	$(CXX) $(CPPFLAGS) -g -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build interp, field dump interpolator.
INTOBJ = interp.o feml.o mesh.o element.o auxfield.o geometry.o

interp: $(INTOBJ)
	$(CXX) $(CPPFLAGS) -g -o $@ $(INTOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build integral, field dump scalar integral evaluator.
INTOBJ = integral.o feml.o mesh.o element.o auxfield.o geometry.o

integral: $(INTOBJ)
	$(CXX) $(CPPFLAGS) -g -o $@ $(INTOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build probe, field dump data extraction.
PRBOBJ = probe.o feml.o mesh.o element.o auxfield.o geometry.o

probe: $(PRBOBJ)
	$(CXX) $(CPPFLAGS) -g -o $@ $(PRBOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build rstress, Reynolds stress calculator.

rstress: rstress.c
	$(CC) $(CPPFLAGS) -o $@ $@.c $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build avgdump, field average updater.

avgdump: avgdump.c
	$(CC) $(CPPFLAGS) -o $@ $@.c $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build rectmesh, rectangular mesh generator;

rectmesh: rectmesh.C
	$(CXX) $(CPPFLAGS) -o $@ $@.C -g

# ----------------------------------------------------------------------------
#
$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Special machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~ *.ti
	rm -rf ILDUMPS
	rm -rf ii_files
	rm -rf ptrepository
