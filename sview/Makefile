##############################################################################
# Makefile for sview, an isosurface visualization package for OpenGL.
#
# $Id$
#

ifdef DEBUG
  DEFINES  = -DDEBUG
else
  DEFINES  =
endif

X11_ROOT  = /usr/X11R6
GLUT_ROOT = /usr
SGI_GLUT  = -lglut -lGLU -lGL -lXmu -lXi -lXext -lX11

CXX      = g++
LD       = g++
CPPFLAGS = -I. -I$(GLUT_ROOT)/include
CXXFLAGS = $(DEFINES) -gstabs+
LDFLAGS  = -L$(X11_ROOT)/lib -L$(GLUT_ROOT)/lib $(SGI_GLUT) -lm

HDRS     = Sview.h
SOBJ     = main.o graphics.o keycom.o semIO.o

sview_3:	$(SOBJ)
	$(LD) -o $@ $(SOBJ) $(LDFLAGS)

clean:
	rm -f *.o

$(SOBJ):	$(HDRS)
