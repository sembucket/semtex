##############################################################################
# Makefile for spectral element solver utilities.
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

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x
CC       = cc
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -fullwarn -mips2 -g
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -lfpe -mips2 -L$(SEM)/lib/$(ARCH) -lfem -lalp \
	   $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CC       = cc
CXX	 = CC
LD	 = CC
CXXFLAGS =  $(DEFINES) -mips4 -64 -Ofast=IP30 -r10000
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -mips4 -64 -Ofast=IP30 -r10000 -L$(SEM)/lib/$(ARCH) -lfem -lalp \
           $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CC       = cc
CXX	 = cxx
LD	 = cxx
CXXFLAGS =  $(DEFINES) -migrate -tune host -g
F77LIBS  = -ldxml
LDFLAGS  = $(CXXFLAGS) -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CC	 = gcc
CXX	 = g++
LD	 = g++
CXXFLAGS = $(DEFINES) -O2
F77LIBS  = -llapack -lblas -lg2c
LDFLAGS  = -O2 -L$(SEM)/lib/$(ARCH) -L/usr/local/lib	\
	   -lfem -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),UNIX_System_V)    # -- Fujitsu VPP300.
LIBDIR   = -L$(SEM)/lib/$(ARCH) 
LIBDIR  += -L/opt/blas/blas_vpp300/lib -L/opt/LAPACK/lib -L/usr/uxplib
LIBS     = -lfem -lalp
LIBS    += -llapack -lblasvpp  -lm -lfj90fv -lfj90 -lfj90f  -ljsp -lvfl
CC       = cc
CXX      = /usr/uxp/C++/bin/CC
LD       = /usr/uxp/C++/bin/CC
#OPTIM    = -Wv,-m3,-Ps 
OPTIM    = -g
CXXFLAGS = $(DEFINES)  $(OPTIM)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS)
endif

ifeq ($(ARCH),SUPER-UX)         # -- NEC SX-4.
LIBDIR   = -L$(SEM)/lib/$(ARCH) 
LIBDIR  += 
LIBS     = -lfem -llapack_64 -lalp -lblas_64 -lfft_64	\
           -lu77sxe -lp77sx -lv77sxe -lf77sxe -li77sx -lm
CC       = cc
CXX      = CC
LD       = CC
#OPTIM    = -g -hacct
CFLAGS   = $(DEFINES) -pvctl,fullmsg -hacct -Onooverlap        \
           -hfloat0 -hnostkchk
CXXFLAGS = $(CFLAGS) $(OPTIM) +a1 -ispeed
CXXFLAGS = $(DEFINES)  $(OPTIM)
LDFLAGS  = $(OPTIM) $(CPPFLAGS) $(LIBDIR) $(LIBS) -hacct
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
# Build enumerate utility.
ENUOBJ = feml.o mesh.o newmerate.o

newmerate:	headers $(ENUOBJ)
	$(LD) -o $@ $(ENUOBJ) $(LDFLAGS) -L. -lmetis

# ----------------------------------------------------------------------------
# Build calc utility.

calc: calc.o
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build addvort utility.
ADDOBJ = addvort.o auxfield.o BCmgr.o boundary.o domain.o	\
	 element.o feml.o field.o geometry.o  integration.o mesh.o	\
	 misc.o	pressure.o

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
ADFOBJ = addfield.o auxfield.o BCmgr.o boundary.o domain.o	\
	 element.o feml.o field.o geometry.o  integration.o mesh.o	\
	 misc.o	pressure.o

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
COMOBJ = compare.o feml.o mesh.o element.o auxfield.o	\
	 geometry.o pressure.o integration.o

compare: $(COMOBJ)
	$(LD) -o $@ $(COMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build noiz, velocity field perturber.

noiz: noiz.c
	$(CC) $(CPPFLAGS) -o $@ $@.c $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build project, field dump projection to higher/lower orders.

project: project.C
	$(CXX) $(CPPFLAGS) -g -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build transform, which does forward/inverse Fourier/Legendre transforms.
TRNOBJ = transform.o auxfield.o element.o geometry.o mesh.o feml.o

transform: transform.C
	$(CXX) $(CPPFLAGS) -g -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build interp, field dump interpolator.
INTOBJ = interp.o feml.o mesh.o element.o auxfield.o geometry.o

interp: $(INTOBJ)
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
	$(CXX) -I$(SEM)/include -c -hnostkchk +a1 -hacct -hnovector -p \
	$(VPATH)/mesh.C
BCmgr.o: $(VPATH)/BCmgr.C
	$(CXX) -I$(SEM)/include -c -hnostkchk +a1 -hacct -hnovector -p \
	$(VPATH)/BCmgr.C
endif

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
	rm -rf ptrepository
