##############################################################################
# Makefile for sview, an isosurface visualisation package for OpenGL.
#
# $Id$
#

ARCH  = $(shell uname -s)

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
endif

X11_ROOT  = /usr/X11R6
GLUT_ROOT = /usr
SGI_GLUT  = -lglut -lGLU -lGL -lXmu -lXi -lXext -lX11

# -- Select compilation system.

ifeq ($(ARCH),IRIX64)		# -- SGI 64-bit systems.
CXX	 = CC
LD	 = CC
CPPFLAGS = -I. -I$(GLUT_ROOT)/include
CXXFLAGS = $(DEFINES) -64 -mips4 -g
LDFLAGS  = -64 -mips4 -L$(X11_ROOT)/lib -L/usr/lib64 \
           -L$(GLUT_ROOT)/lib $(SGI_GLUT) -lm
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CXX      = g++
LD       = g++
CPPFLAGS = -I. -I$(GLUT_ROOT)/include
CXXFLAGS = $(DEFINES) -O3
LDFLAGS  = -L$(X11_ROOT)/lib -L$(GLUT_ROOT)/lib $(SGI_GLUT) -lm
endif

# -- Remainder is machine-independent.

HDRS     = Sview.h
SOBJ     = main.o graphics.o isowire.o keycom.o semIO.o

sview:	$(SOBJ)
	$(LD) -o $@ $(SOBJ) $(LDFLAGS)

clean:
	rm -f *.o

$(SOBJ):	$(HDRS)
