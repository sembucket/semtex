##############################################################################
# Makefile for arnoldi eigensolver.
#
# $Id$
##############################################################################

SEM      = ../semtex-5.2
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
FFLAGS   = $(CXXFLAGS)
CFLAGS   = $(CXXFLAGS)
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -lfpe -mips2 -L$(SEM)/lib/$(ARCH) -lfem -lalp \
	   $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CC       = cc
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -g -mips4 -64
FFLAGS   = $(CXXFLAGS)
CFLAGS   = $(CXXFLAGS)
F77LIBS	 = -lcomplib.sgimath -lftn
LDFLAGS  = -mips4 -64 -L$(SEM)/lib/$(ARCH) -lalp $(F77LIBS) -lmalloc -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CC       = cc
CXX	 = cxx
LD	 = cxx
CXXFLAGS = $(DEFINES) -migrate -tune host -g
FFLAGS   = $(CXXFLAGS)
CFLAGS   = $(CXXFLAGS)
F77LIBS  = -ldxml
LDFLAGS  = $(CXXFLAGS) -L$(SEM)/lib/$(ARCH) -lalp $(F77LIBS) -lm
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CC	 = gcc
CXX	 = g++
LD	 = g++
CXXFLAGS = $(DEFINES) -gstabs+
FFLAGS   = $(CXXFLAGS)
CFLAGS   = $(CXXFLAGS)
F77LIBS  = -llapack -lblas -lf2c
LDFLAGS  = -gstabs+ -L$(SEM)/lib/$(ARCH) -L/usr/local/lib	\
	   -lalp $(F77LIBS) -lm
endif

# ----------------------------------------------------------------------------
# Build arnoldi solver:

OBJ = iohb.o dspmvc.o arnoldi.o

arnoldi:	$(OBJ)
	$(LD) -o $@ $(OBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Admin stuff.

clean:
	rm -f *.o *~
	rm -rf ILDUMPS ii_files ptrepository
