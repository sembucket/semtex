##############################################################################
# Makefile for qmesh 2D mesh generation.
#
# $Id$
##############################################################################

ARCH     = $(shell uname -s)
CPPFLAGS = -I. -I$(HOME)/include -I$(HOME)/include/sm
DEFINES  = 

# -- Set DEBUG=1 on command line if you want debugging output.

ifdef DEBUG
  DEFINES += -DDEBUG
endif

# -- Set GRAPHICS=1 on command line if you want to compile in SM calls.

ifdef GRAPHICS
  DEFINES += -DGRAPHICS
endif

# -- Set PROMPT=1 on command line if you want loop-by-loop prompting

ifdef PROMPT
  DEFINES += -DPROMPT
endif

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x

SMLIBS    = -lplotsub -ldevices -lutils
X11LIBS   = -lX11_s -lgl_s

CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -mips2 -fullwarn -smart -g
LDFLAGS  = -L$(HOME)/lib -L/usr/local/lib/sm -mips2 -lm
endif


ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -O3
LDFLAGS  = $(CXXFLAGS) -lm
endif


ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CXX	 = cxx
LD	 = cxx
CXXFLAGS = $(DEFINES) -migrate -tune host -O2
LDFLAGS  = $(CXXFLAGS) -lm
endif


# ----------------------------------------------------------------------------
#
SRC = main loop node point graphics
OBJ = $(addsuffix .o,$(SRC))
HDR = qmesh.h

# ----------------------------------------------------------------------------
#

qmesh: $(HDR) $(OBJ)
	$(LD) -o qmesh $(OBJ) $(CPPFLAGS) $(LDFLAGS) $(SMLIBS) $(X11LIBS)

# ----------------------------------------------------------------------------
#
$(OBJ): $(HDR)

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
