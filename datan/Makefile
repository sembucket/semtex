OPTIONS =
CC      = cc
CFLAGS  = $(OPTIONS)
CPPFLAGS= $(DEFINES) -I$(HOME)/include
INCLUDES= $(HOME)/include/nrutil.h $(HOME)/include/FFTutil.h 
LD      = $(CC)
LDFLAGS = $(CFLAGS)
LIBS    = -L$(HOME)/lib -lFFTutil -lnrutil -lm

MEM     = memspec
SRCS    = main.c io.c memcof.c evlmem.c
OBJS    = main.o io.o memcof.o evlmem.o

acor: acor.c $(INCLUDES)
	$(LD) $(LDFLAGS) $(CPPFLAGS) -o acor acor.c $(LIBS)

chop: chop.c
	$(LD) $(LDFLAGS) -o chop chop.c

pdf: pdf.c $(HOME)/include/nrutil.h
	$(LD) $(LDFLAGS) $(CPPFLAGS) -o pdf pdf.c $(LIBS)

slit: slit.c $(INCLUDES)
	$(LD) $(LDFLAGS) $(CPPFLAGS) -o slit slit.c $(LIBS)

testpow: testpow.c $(INCLUDES)
	$(LD) $(LDFLAGS) $(CPPFLAGS) -o testpow testpow.c $(LIBS)

clean:
	rm -rf *.o *~ core
