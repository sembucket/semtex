##############################################################################
# Makefile for "fem++" finite element solver.
#
# $Id$
##############################################################################

ARCH     = $(shell uname -s)
PROG     = fem++

OPTIONS  = -O2
DEFINES  =

CXX      = CC
LD       = CC

CXXFLAGS = $(OPTIONS) $(DEFINES)
CPPFLAGS = -I../include
F77LIBS  = -lcomplib.sgimath -lisam
LDFLAGS  = -L../lib -lfem -lalp $(F77LIBS) -lm -lfpe

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)
override CXXFLAGS += -mips2 -fullwarn
override LDFLAGS  += -mips2
endif

ifeq ($(ARCH),IRIX64)

endif

# ----------------------------------------------------------------------------


FILES    = drive domain field element BC boundary pressure elliptic NS misc
OBJS     = $(addsuffix .o,$(FILES))
HEADER   = Fem.h



$(PROG): $(HEADER) $(OBJS)
	$(LD) -o $(PROG).tmp $(OBJS) $(CPPFLAGS) $(LDFLAGS)
	mv $(PROG).tmp $(PROG)

$(OBJS): $(HEADER)

clean:
	rm -f *.o *~







