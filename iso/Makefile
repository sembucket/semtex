##############################################################################
# This directory contains code for simulation of isotropic turbulence by a
# Fourier--Galerkin pseudo-spectral method.   See more extensive comments
# in file drive.c
#
# Hugh Blackburn
# Department of Mechanical Engineering
# Monash University
# Clayton 3168
# Australia
# hmb@artemis.eng.monash.edu.au
#
# $Id$
##############################################################################
# Compilation options and flags.

ISO      = .
VPATH    = $(ISO)/src

ARCH     = $(shell uname -s)
PROG     = iso

DEFINES  =

CPPFLAGS = -I.

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x
CC	 = cc
LD	 = cc
CFLAGS   = -mips2 -fullwarn -g
LDFLAGS  = $(CFLAGS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CC       = cc
LD       = cc
CFLAGS   = -O3
LDFLAGS  = $(CFLAGS) -lmalloc -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1 Alpha
CC       = cc
LD	 = cc
CFLAGS   = -O2 -migrate -tune host
LDFLAGS  = $(CFLAGS) -lm
endif

HEADER    = iso.h

# ----------------------------------------------------------------------------
# Make iso application.
ISOOBJ    = allocate.o io.o nonlinear.o integrate.o energy.o FFT.o taylor.o drive.o

$(PROG): $(ISOOBJ) $(HEADER)
	$(LD) $(ISOOBJ) -o $@ $(LDFLAGS)

# ----------------------------------------------------------------------------
# Make a top-hat spectrum initial conditions file.
IC1OBJ = allocate.o io.o normalize.o energy.o ran2PI.o tophat.o FFT.o makeIC1.o

makeIC1: $(IC1OBJ) $(HEADER)
	$(LD) $(IC1OBJ) -o $@ $(LDFLAGS)

# ----------------------------------------------------------------------------
# Make ICs for 3D Taylor-Green vortex.
IC2OBJ = allocate.o io.o FFT.o taylor.o makeIC2.o

makeIC2: $(IC2OBJ) $(HEADER)
	$(LD) $(IC2OBJ) -o $@ $(LDFLAGS)

# ----------------------------------------------------------------------------
# Make ICs for 3D Taylor-Green vortex.
IC3OBJ = allocate.o io.o FFT.o taylor.o makeIC3.o
IC3CHK = allocate.o io.o energy.o FFT.o taylor.o taylor2D_chk.o

makeIC3: $(IC3OBJ) $(IC3CHK) $(HEADER)
	$(LD) $(IC3OBJ) -o $@ $(LDFLAGS)
	$(LD) $(IC3CHK) -o taylor2D_chk $(LDFLAGS)

# ----------------------------------------------------------------------------
# Test the FFT routines.
TESTOBJ = allocate.o io.o FFT.o testDFT.o

testDFT: $(TESTOBJ) $(HEADER)
	$(CC) -o $@ $(TESTOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Divergence-checking diagnostic.
CHKOBJ = allocate.o io.o derivative.o FFT.o truncation.o div_chk.o

div_chk: $(CHKOBJ) $(HEADER)
	$(CC) -o $@ $(CHKOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Convert field file to tecplot format, with vorticity too.
TECOBJ = allocate.o io.o derivative.o truncation.o FFT.o iso2tec.o

iso2tec: $(TECOBJ) $(HEADER)
	$(LD) $(CFLAGS) -o $@ $(TECOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Clean up.
clean:
	rm -f *.o *.u *~



