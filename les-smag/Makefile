##############################################################################
# Makefile for finite element solver.
#
# $Id$
##############################################################################

FEM      = ..
VPATH    = $(FEM)/src

ARCH     = $(shell uname -s)

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
endif

CPPFLAGS = -I. -I$(FEM)/include

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -mips2 -fullwarn -smart -g
F77LIBS	 = -lcomplib.sgimath -lisam
LDFLAGS  = -lfpe -mips2 -L$(FEM)/lib/$(ARCH) -lfem -lalp \
	   $(F77LIBS) -lmalloc -lm
endif


ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CXX	 = CC
LD	 = CC
CXXFLAGS =  $(DEFINES) -O3
F77LIBS	 = -lcomplib.sgimath -lisam
LDFLAGS  = -L$(FEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lmalloc -lm
endif


ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CXX	 = cxx
LD	 = cxx
CXXFLAGS =  $(DEFINES) -migrate -tune host -O2
F77LIBS  = -ldxml
LDFLAGS  = $(CXXFLAGS) -L$(FEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif


# ----------------------------------------------------------------------------
#
FEMFILES  = BC boundary domain element family field mesh misc pressure
FEMOBJ    = $(addsuffix .o,$(FEMFILES))
FEMHDR    = Fem.h

# ----------------------------------------------------------------------------
# Default build for Navier-Stokes solver.
NS_OBJ = drive.o NS.o

ns:	$(FEMOBJ) $(NS_OBJ)
	$(LD) -o $@ $(NS_OBJ) $(FEMOBJ) $(LDFLAGS)
# ----------------------------------------------------------------------------
#
$(FEMOBJ): $(FEMHDR)

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
