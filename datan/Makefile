OPTIONS = -g -64
CC      = cc
CFLAGS  = $(OPTIONS)
CPPFLAGS= $(DEFINES)
LDFLAGS = -L/home/hmb/lib/IRIX64 -lrecipesDP
INCLUDE = FFTutil.h

asp: FFTutil.o asp.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o asp FFTutil.o asp.o -lm

acor: FFTutil.o acor.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o acor FFTutil.o acor.o -lm

correlate: FFTutil.o correlate.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ FFTutil.o $@.o $(LDFLAGS) -lm

chop: chop.c
	$(CC) $(CFLAGS) -o chop chop.c

slit: FFTutil.o slit.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o slit FFTutil.o slit.o -lm

pdf: FFTutil.o pdf.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o pdf FFTutil.o pdf.o -lm

memspec: FFTutil.o memspec.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o memspec FFTutil.o memspec.o -lm

dft_chk: FFTutil.o dft_chk.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o dft_chk FFTutil.o dft_chk.o -lm

pcor: pcor.c
	$(CC) $(CFLAGS) -o pcor pcor.c
clean:
	rm -rf *.o *~ core
