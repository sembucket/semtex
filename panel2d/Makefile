#############################################################################
# Makefile for panel-method codes.
# RCS Information:
# $Id$
#############################################################################

OPTIONS  = -g

CFLAGS   = $(OPTIONS)
CPPFLAGS = -I$(HOME)/include
LDFLAGS  = -lm


OBJ       = panel2d.o lu.o singular.o
HEADER    = panel.h

panel2d: $(OBJ) $(HEADER)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

OBJ: $(HEADER)

clean:
	-rm -f *.o a.out *~ core




