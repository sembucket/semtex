##############################################################################
# This directory contains code for simulation of isotropic turbulence by a
# Fourier Galerkin pseudo-spectral method.   See more extensive comments in 
# the driver file drive.c
# $Id$


##############################################################################
# Compile for SGI 4D and INDIGO running IRIX.R4
.SUFFIXES:
.SUFFIXES: .c .o
DEBUGF    = -g -fullwarn -woff 269-270,287 -float
DEBUGN    = iso.db
FINALF    = -O3 -float
FINALN    = iso
LIBS      = -lm
CC        = cc

HEADER    = globals.h

.c.o:; cc -c -o $@ $(DEBUGF) $*.c

##############################################################################
# Make iso.db with debugging compilation. (DEFAULT).
ISOOBJ    = allocate.o io.o convolve.o integrate.o nonlin.o rc3DFT.o\
            rcZFFT.o cZFFT.o cYFFT.o cXFFT.o pcYFFT.o pcXFFT.o\
            maketables.o drive.o energy.o

iso.db: $(ISOOBJ) $(HEADER)
	$(CC) $(DEBUGF) -o $(DEBUGN) $(ISOOBJ) $(LIBS)


##############################################################################
# Make iso with -O3 optimization, after debugging.
ISOSRC    = allocate.c io.c convolve.c integrate.c nonlin.c rc3DFT.c\
            rcZFFT.c cZFFT.c cYFFT.c cXFFT.c pcYFFT.c pcXFFT.c\
            maketables.c drive.c energy.c

iso: 	$(ISOSRC) $(HEADER)
	$(CC) $(FINALF) -s -o $(FINALN) $(ISOSRC) $(LIBS)


##############################################################################
# Make a top-hat spectrum initial conditions file.
IC1OBJ = allocate.o io.o normalize.o ran2PI.o tophat.o rc3DFT.o\
         rcZFFT.o cZFFT.o cYFFT.o cXFFT.o pcYFFT.o pcXFFT.o\
	 makeIC1.o energy.o maketables.o

makeIC: $(IC1OBJ) $(HEADER)
	$(CC) $(DEBUGF) -o makeIC $(IC1OBJ) $(LIBS)


##############################################################################
# Test the FFT routines.
TESTOBJ = allocate.o io.o testDFT.o maketables.o rc3DFT.o\
          rcZFFT.o cZFFT.o cYFFT.o cXFFT.o pcYFFT.o pcXFFT.o

testDFT:	$(TESTOBJ) $(HEADER)
	$(CC) $(DEBUGF) -o $@ $(TESTOBJ) $(LIBS)

##############################################################################
# Divergence-checking diagnostic.
CHKOBJ = allocate.o io.o maketables.o deriv.o div_chk.o rc3DFT.o\
         rcZFFT.o cZFFT.o cYFFT.o cXFFT.o pcYFFT.o pcXFFT.o

div_chk: $(CHKOBJ) $(HEADER)
	$(CC) $(DEBUGF) -o $@ $(CHKOBJ) $(LIBS)

##############################################################################
# Clean up.
clean:
	rm *.o *.u
