##############################################################################
# Makefile for "fem++" finite element solver.
#
# $Id$
##############################################################################

# ----------------------------------------------------------------------------
# If you are compiling from a remote directory, set FEM to be the
# name of the directory where the fem++ source code and libraries are 
# installed

FEM      = .
VPATH    = $(FEM)/src

ARCH     = $(shell uname -s)
PROG     = fem++

OPTIONS  = -g
DEFINES  =

CXX      = CC
LD       = CC

CXXFLAGS = $(OPTIONS) $(DEFINES)
CPPFLAGS = -I$(FEM)/include
F77LIBS  = -lcomplib.sgimath -lisam
LDFLAGS  = -L$(FEM)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)
override CXXFLAGS += -mips2 -fullwarn -smart
override LDFLAGS  += -lfpe -mips2
endif

ifeq ($(ARCH),IRIX64)

endif

# ----------------------------------------------------------------------------
# Object code and header files required by the program.
#
# If you have any extra source code to compile, add it to the list below,
# preferably as the variable SPECIAL.  The source code postfix (assumed .C) 
# and object name (assumed .o) are generated automatically.

FILES   = drive domain field element BC boundary pressure elliptic NS misc
SPECIAL =
OBJS    = $(addsuffix .o,$(FILES))
HEADER  = Fem.h


$(PROG): $(HEADER) $(OBJS)
	$(LD) -o $(PROG).tmp $(OBJS) $(CPPFLAGS) $(LDFLAGS)
	mv $(PROG).tmp $(PROG)

$(OBJS): $(HEADER)

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
		util/*.C util/*.c util/*.h examples
	compress -f dist.tar

# ----------------------------------------------------------------------------
# Make full source code distribution with library sources.

fulldist:
	tar cvf fulldist.tar Makefile include/*.h src/*.C src/*.h lib	\
		util/*.C util/*.c util/*.h examples alplib femlib
	compress -f fulldist.tar

# ----------------------------------------------------------------------------
# Clean up.

clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
