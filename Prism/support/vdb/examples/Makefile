#
# Makefile for the VDB examples
#

include $(PRISM)/system/defaults

PROGS    = hello test dot sum mxv
CPPFLAGS = -I.. -I$(PRISM)/include -g
LIBS     = ../libvdb.a -L$(PRISM)/lib $(CUBIT) $(VECLIB) \
	$(BLAS) $(LAPACK) $(COMM) $(SYSTEM)

# ------------------------------------------------------------------

all: $(PROGS)

hello: hello.c
	$(CC) $(CPPFLAGS) -o hello hello.c $(LIBS)

test: test.c
	$(CC) $(CPPFLAGS) -o test test.c $(LIBS)

dot: dot.c
	$(CC) $(CPPFLAGS) -o dot dot.c $(LIBS)

sum: sum.c
	$(CC) $(CPPFLAGS) -o sum sum.c $(LIBS)

mxv: mxv.c
	$(CC) $(CPPFLAGS) -o mxv mxv.c $(LIBS)

# ------------------------------------------------------------------

tidy:
	$(RM) *~ \#*\#
clean: tidy
	$(RM) *.o core* $(PROGS)
empty: clean
