##############################################################################
# Makefile for "fem++" finite element solver.
#
# $Id$
##############################################################################

# ----------------------------------------------------------------------------
# If you are compiling from a remote directory, set FEM to be the
# name of the directory where the fem++ source code and library
# directories are installed.

FEM      = .
VPATH    = $(FEM)/src

ARCH     = $(shell uname -s)
PROG     = fem++

OPTIONS  = -g
DEFINES  =

CPPFLAGS = -I. -I$(FEM)/include

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x
CXX	 = CC
LD	 = CC
CXXFLAGS = -mips2 -fullwarn -smart
F77LIBS	 = -lcomplib.sgimath -lisam
LDFLAGS  = -lfpe -mips2 -L$(FEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif


ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CXX	 = CC
LD	 = CC
CXXFLAGS = -O3
F77LIBS	 = -lcomplib.sgimath -lisam
LDFLAGS  = -L$(FEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif


ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CXX	 = cxx
LD	 = cxx
CXXFLAGS = -O2 -migrate -tune host
F77LIBS  = -ldxml
LDFLAGS  = $(CXXFLAGS) -L$(FEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
endif


# ----------------------------------------------------------------------------
# Default build for fem++.
#
# If you have any extra source code to compile, add it to the list below,
# preferably as the variable SPECIAL.  The source code postfix (assumed .C) 
# and object name (assumed .o) are generated automatically.

FEMFILES  = drive     domain   field   element  BC     boundary   pressure \
            elliptic  NS       misc    mesh
SPECIAL   =
FEMOBJ    = $(addsuffix .o,$(FEMFILES))
FEMHDR    = Fem.h


$(PROG): $(FEMHDR) $(FEMOBJ)
	$(LD) -o $(PROG).tmp $(FEMOBJ) $(CPPFLAGS) $(LDFLAGS)
	mv $(PROG).tmp $(PROG)

$(FEMOBJ): $(FEMHDR)

# ----------------------------------------------------------------------------
# Build meshpr utility.

meshpr:	meshpr.C meshpr.o field.o mesh.o element.o misc.o boundary.o BC.o \
	pressure.o
	$(LD) $(CXXFLAGS) meshpr.o field.o element.o mesh.o misc.o \
		boundary.o BC.o pressure.o -o $@ $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build calc utility.

calc:	calc.C
	$(LD) $(CFLAGS) -o $@ $@.C $(CPPFLAGS) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Run this to compile libraries required by program.

libs:
	cd $(FEM)/alplib;	\
	gmake install;		\
	gmake;			\
	gmake install

	cd $(FEM)/femlib;	\
	gmake install;		\
	gmake;			\
	gmake install

# ----------------------------------------------------------------------------
# Make base distribution with source code and compiled libraries.

dist:
	tar cvf dist.tar Makefile include/*.h src/*.C src/*.h lib	\
		util/*.C util/*.c util/*.h examples bin
	compress -f dist.tar

# ----------------------------------------------------------------------------
# Make full source code distribution with library sources.

fulldist:
	tar cvf fulldist.tar Makefile include/*.h src/*.C src/*.h lib	\
		util/*.C util/*.c util/*.h examples alplib femlib bin
	compress -f fulldist.tar

# ----------------------------------------------------------------------------
# Clean up.

clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
