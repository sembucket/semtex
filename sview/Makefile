##############################################################################
# Makefile for sview, an isosurface visualisation tool for OpenGL.
#
# Copyright (C) 1999 Hugh Blackburn
#
# Set DEBUG=1  when making to make a debuggable executable.
# Set PYTHON=1 when making to build bindings for python.
#
# $Id$
#

ARCH  = $(shell uname -s)

# -- These may need editing: they specify where to find X & GLUT libraries:

X11_ROOT  = /usr/X11R6
GLUT_ROOT = /usr

# -- If compiling Python bindings tell the system where the headers are:

PYTHON_HDRS = /usr/local/include/python1.3

# -- Switch on any conditional compilation for debugging:

ifdef DEBUG
  DEFINES = -DDEBUG
else
  DEFINES =
endif

# -- Select compilation system.

ifeq ($(ARCH),IRIX64)		# -- SGI 64-bit systems.

CXX = CC
LD  = CC

GLUT_LIBS = -ltiff -lglut -lGLU -lGL -lX11 -lXext -lXmu -lXt -lSM -lICE

CPPFLAGS = -I. -I$(GLUT_ROOT)/include
ifdef PYTHON
  CPPFLAGS += -I$(PYTHON_HDRS)
endif

CXXFLAGS = $(DEFINES)
ifdef DEBUG
  CXXFLAGS += -g -32
else
  CXXFLAGS += -O2 -32
endif

LDFLAGS = -32 -L$(X11_ROOT)/lib -L$(GLUT_ROOT)/lib $(GLUT_LIBS) -lm
endif


ifeq ($(ARCH),Linux)		# -- GNU/Linux system.

CXX = g++
LD  = g++

GLUT_LIBS = -ltiff -ljpeg -lz -lglut -lGLU -lGL	\
            -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

CPPFLAGS = -I. -I$(GLUT_ROOT)/include
ifdef PYTHON
  CPPFLAGS += -I$(PYTHON_HDRS)
endif

CXXFLAGS = $(DEFINES)
ifdef DEBUG
  CXXFLAGS += -gstabs+
else
  CXXFLAGS += -O3
endif

LDFLAGS  = -L$(X11_ROOT)/lib -L$(GLUT_ROOT)/lib $(GLUT_LIBS) -lm
endif


# -- Remainder is machine-independent.

HDRS = ./Sview.h
SOBJ = main.o graphics.o image.o isowire.o keycom.o semIO.o

DEPEND = $(SOBJ)
ifdef PYTHON
  DEPEND += Sview.so
endif

sview: $(DEPEND)
	$(LD) $(SOBJ) $(LDFLAGS) -o $@

$(SOBJ): $(HDRS)

Sview.so: bindings.o graphics.o isowire.o keycom.o semIO.o
	$(LD) $(LDFLAGS) \
	-shared bindings.o graphics.o isowire.o keycom.o semIO.o -o $@
tidy:
	$(RM) *~ \#*\#
clean: tidy
	$(RM) -r *.o *.so core ii_files so_locations
empty: clean
	$(RM) sview
