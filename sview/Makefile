##############################################################################
# Makefile for sview, an isosurface visualisation tool for OpenGL.
#
# Copyright (c) 1999 <--> $Date$, Hugh Blackburn
#
# Set DEBUG=1  when making to make a debuggable executable.
# Set PYTHON=1 when making to build bindings for python.
#
# $Id$
#

ARCH  = $(shell uname -s)
MACH  = $(shell uname -m)

# -- These may need editing: they specify where to find X & GLUT libraries:

X11_ROOT  = /usr/X11R6
ifeq ($(MACH),x86_64)
  GLUT_ROOT = /usr/local
else
  GLUT_ROOT = /usr
endif

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
  CXXFLAGS += -g -n32
else
  CXXFLAGS += -O2 -n32
endif
LDFLAGS = -n32 -L$(X11_ROOT)/lib -L$(GLUT_ROOT)/lib32 $(GLUT_LIBS) -lm
endif


ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CXX = g++
LD  = g++
GLUT_LIBS = -ltiff -ljpeg -lz -lglut -lGLU -lGL	\
            -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE
CPPFLAGS = -I.
ifdef PYTHON
  CPPFLAGS += -I$(PYTHON_HDRS)
endif
CXXFLAGS = $(DEFINES)
ifdef DEBUG
  CXXFLAGS += -gstabs+
else
  CXXFLAGS += -O2
endif
LDFLAGS  = -L$(X11_ROOT)/lib -L$(GLUT_ROOT)/lib $(GLUT_LIBS) -lm
endif

ifeq ($(MACH),x86_64)		# -- GNU/Linux/Opteron system.
CXX = g++
LD  = g++
GLUT_LIBS = -ltiff -ljpeg -lz /usr/lib64/libglut.a -lGLU -lGL 		\
		-lX11 -lXext -lXmu -lXt -lXi -lSM -lICE
CPPFLAGS = -I.
ifdef PYTHON
  CPPFLAGS += -I$(PYTHON_HDRS)
endif
CXXFLAGS = $(DEFINES)
ifdef DEBUG
  CXXFLAGS += -gstabs+
else
  CXXFLAGS += -O2
endif
LDFLAGS  = -L$(X11_ROOT)/lib64 -L$(GLUT_ROOT)/lib $(GLUT_LIBS) -lm
endif


ifeq ($(ARCH),Darwin)		# -- Mac OSX (10.3, 10.4).
CXX = g++
LD  = g++
CPPFLAGS = -I. -I/sw/include
ifdef PYTHON
  CPPFLAGS += -I$(PYTHON_HDRS)
endif
CXXFLAGS = $(DEFINES)
ifdef DEBUG
  CXXFLAGS += -gstabs+
else
  CXXFLAGS += -O2
endif
# -- Using the OSX GLUT framework in OSX 10.3 gives problems with command-line
#    filename resolution. Instead, we use the fink-supplied libglut:
#LDFLAGS = -framework OpenGL -L/sw/lib -ltiff -lglut
# -- But in OSX 10.4 it works fine:
LDFLAGS = -framework OpenGL -framework GLUT -L/sw/lib -ltiff
endif


# -- Remainder is machine-independent.

HDRS = ./sview.h
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
