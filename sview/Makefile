##############################################################################
# Makefile for SVIEW_2, an isosurface visualization package for SGI machines.
#
# $Id$
#
CC	= cc
LD	= cc
INCDIR  = -I. -I/usr/local/include
LFLAGS 	= -O2 
CFLAGS	= -O2 -DPI $(INCDIR)
LIBS    = -lgl -lm -lsun -limage -lbsd -lmalloc -lgutil

SOBJ    = main_vis.o display_surface.o flow_io.o m_cube.o	\
	  make_mesh.o misc_make.o misc_view.o

sview_2:	$(SOBJ)
	$(LD) $(CFLAGS) $(LFLAGS) -o $@ $(SOBJ) $(LIBS)

clean:
	rm -f *.o
