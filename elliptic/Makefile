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
CXXFLAGS = $(DEFINES) -fullwarn -smart -g
F77LIBS	 = -lcomplib.sgimath -lisam
LDFLAGS  = -lfpe -L$(SEM)/lib/$(ARCH) -lfem -lalp \
	   $(F77LIBS) -lmalloc -lmx
endif

ifeq ($(ARCH),IRIX64)		# -- SGI 64-bit systems.
CXX	 = CC
LD	 = CC
CXXFLAGS =  $(DEFINES) -g
F77LIBS	 = -lcomplib.sgimath -lisam
LDFLAGS  = -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CXX	 = cxx
LD	 = cxx
CXXFLAGS =  $(DEFINES) -migrate -tune host -g -w0
F77LIBS  = -ldxml
LDFLAGS  = $(CXXFLAGS) -L$(SEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CXX	 = g++
LD	 = g++
CXXFLAGS = $(DEFINES) -gstabs+
F77LIBS  = -llapack -lblas -lf2c
LDFLAGS  = -gstabs+ -L$(SEM)/lib/$(ARCH) -L/usr/local/lib	\
	   -lfem -lalp $(F77LIBS)
endif

# ----------------------------------------------------------------------------
#

SEMFILES  = auxfield BCmgr boundary condition domain element	\
	    feml field integration matrix mesh misc pressure
SEMOBJ    = $(addsuffix .o,$(SEMFILES))
SEMHDR    = Sem.h


# ----------------------------------------------------------------------------
#
ELLOBJ = drive.o helmholtz.o

elliptic:	sem $(SEMHDR) $(SEMOBJ) $(ELLOBJ)
	cd ../src; $(MAKE)
	$(LD) -o elliptic $(ELLOBJ) $(SEMOBJ) $(CPPFLAGS) $(LDFLAGS)

# ----------------------------------------------------------------------------
#
$(SEMOBJ): $(SEMHDR)

sem:
	cd $(VPATH); make install

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
