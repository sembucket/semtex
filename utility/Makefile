#############################################################################
# Makefile for finite-element utility routines.
#
# RCS Information:
# $Id$
#############################################################################

CC         = cc
F77        = f77
OPTIONS    = -g
CFLAGS     = $(OPTIONS)
INCDIR     = -I/usr/people/hmb/src/Fem/include
LIBDIR     = -L/usr/people/hmb/src/Fem/lib
LIBS       = -lfembase -lm

cavitymesh:	cavitymesh.c
	$(CC) $(CFLAGS) -o cavitymesh cavitymesh.c $(INCDIR) $(LIBDIR) $(LIBS)

connect:	connect.c
	$(CC) $(CFLAGS) -o connect connect.c $(INCDIR) $(LIBDIR) $(LIBS)

clean:
	rm -rf *.o *~

tar:
	cd .. \
	tar cvf util.tar util \
	compress util.tar
