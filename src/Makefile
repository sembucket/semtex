##############################################################################
# Makefile for finite element method solver.
#
# Written for the GNUmake utility by Stallman & McGrath.
##############################################################################

PROG     = fem

OPTIONS  = -g -mips2

CC       = cc
LD       = cc

CFLAGS   = $(OPTIONS)
CPPFLAGS = -I../include
F77LIBS  = -lcomplib.sgimath -lF77 -lI77 -lU77 -lisam
LDFLAGS  = -L../lib -lfem -lalp $(F77LIBS) -lm -lfpe -lmalloc


FILES    = boundary drive edge element field matrix mesh prepost
SRCS     = $(addsuffix .c,$(FILES))
OBJS     = $(addsuffix .o,$(FILES))
HEADER   = $(PROG).h



$(PROG): $(HEADER) $(OBJS)
	$(LD) $(CFLAGS) -o $(PROG).tmp $(OBJS) $(CPPFLAGS) $(LDFLAGS)
	mv $(PROG).tmp $(PROG)

clean:
	rm -f *.o *~


