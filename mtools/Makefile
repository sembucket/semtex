##############################################################################
# Makefile for qmesh 2D mesh generation.
#
# $Id$
##############################################################################

ARCH     = $(shell uname -s)
#CPPFLAGS = -I. -I$(HOME)/develop/semtex/include -I$(HOME)/include/sm
CPPFLAGS = -I. -I$(HOME)/include/sm
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

ifeq ($(ARCH),Darwin)               # -- Mac OSX.
CXX      = g++
LD       = g++
CC       = gcc
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O2
endif
DEFINES += -DmacOSX		# This is for SM graphics
CFLAGS   = $(DEFINES) $(OPT) -ansi
CXXFLAGS = $(DEFINES) $(OPT) -ansi
LIBDIR   += -L/sw/lib -L/usr/X11R6/lib
CPPFLAGS += -I/sw/include
LDFLAGS  = $(OPT) $(LIBDIR) 
CLDFLAGS = $(LDFLAGS)
SMLIBS   = -lsm_plotsub -lsm_devices -lsm_utils -lX11
endif

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x
SMLIBS    = -lplotsub -ldevices -lutils
X11LIBS   = -lX11_s -lgl_s

CXX	 = CC
LD	 = CC
CXXFLAGS = $(DEFINES) -mips2 -fullwarn -smart -g
LDFLAGS  = -L$(HOME)/lib -L/usr/local/lib/sm -mips2 -lm
endif

ifeq ($(ARCH),IRIX64)		# -- IRIX64
CXX	 = CC
LD	 = CC
#CXXFLAGS = $(DEFINES) -n32 -g -LANG:std
CXXFLAGS = $(DEFINES) -32 -g -LANG:std
LDFLAGS  = $(CXXFLAGS)  -L/usr/local/lib/sm -L/usr/X11/lib -lm
SMLIBS   = -lplotsub -ldevices -lutils
X11LIBS  = -lX11
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CXX	 = cxx
LD	 = cxx
CXXFLAGS = $(DEFINES) -migrate -tune host -O2
LDFLAGS  = $(CXXFLAGS) -lm
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CXX	 = g++
LD	 = g++
CXXFLAGS = $(DEFINES) -g
LDFLAGS  = -L$(HOME)/lib -L/usr/local/lib -L/usr/X11/lib
SMLIBS   = -lsm_plotsubSP -lsm_devicesSP -lsm_utilsSP
X11LIBS  = -lX11
endif

# ----------------------------------------------------------------------------
#
HDR = qmesh.h

# ----------------------------------------------------------------------------
# qmesh --- mesh generation utility.

QSRC = qmesh loop node point graphics misc global
QOBJ = $(addsuffix .o,$(QSRC))

qmesh: $(HDR) $(QOBJ)
	$(LD) -o qmesh $(QOBJ) $(LDFLAGS) $(SMLIBS) $(X11LIBS) -lm

$(QOBJ): $(HDR)

# ----------------------------------------------------------------------------
# reflect --- mesh manipulation utility.

RSRC = reflect node point misc global
ROBJ = $(addsuffix .o,$(RSRC))

reflect: $(HDR) $(ROBJ)
	$(LD) -o reflect $(ROBJ) $(CPPFLAGS) $(LDFLAGS)

$(ROBJ): $(HDR)

# ----------------------------------------------------------------------------
# outline --- input generator for qmesh.

OSRC = outline node point misc
OOBJ = $(addsuffix .o,$(OSRC))

outline: $(HDR) $(OOBJ)
	$(LD) -o outline $(OOBJ) $(CPPFLAGS) $(LDFLAGS)

# ----------------------------------------------------------------------------
# smooth --- utility for Laplacian smoothing (hand-built) meshes.

SSRC = smooth node point misc
SOBJ = $(addsuffix .o,$(SSRC))

smooth: $(HDR) $(SOBJ)
	$(LD) -o smooth $(SOBJ) $(CPPFLAGS) $(LDFLAGS) $(SMLIBS) $(X11LIBS)

# ----------------------------------------------------------------------------
# merge --- combine two qmesh input files.

MSRC = merge loop node point misc graphics global
MOBJ = $(addsuffix .o,$(MSRC))

merge: $(HDR) $(MOBJ)
	$(LD) -o merge $(MOBJ) $(CPPFLAGS) $(LDFLAGS) $(SMLIBS) $(X11LIBS) -lm

# ----------------------------------------------------------------------------
# unique --- fix up non-unique nodes produced by digmesh.

UOBJ = unique.o point.o misc.o

unique: $(UOBJ) $(HDR)
	$(LD) -o unique $(UOBJ) $(CPPFLAGS) $(LDFLAGS) -lm


$(OOBJ): $(HDR)

# ----------------------------------------------------------------------------
# Make tarfile of source.

tar:
	tar cvf qmesh.tar *	;	\
	gzip qmesh.tar

# -----------------------------------------------------------------------------
#
clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
