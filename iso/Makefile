##############################################################################
# This directory contains code for simulation of isotropic turbulence by a
# Fourier--Galerkin pseudo-spectral method.   See more extensive comments
# in file drive.c
#
# Hugh Blackburn
# CSIRO
# Building, Construction & Engineering
# PO Box 56
# Highett 3190
# Australia
# hugh.blackburn@dbce.csiro.au
#
# $Id$
##############################################################################
# Compilation options and flags.

ISO      = .
VPATH    = $(ISO)/src

ARCH     = $(shell uname -s)
PROG     = iso

ifdef TG
  DEFINES = -DTG
else
  DEFINES =
endif

ifdef INVISCID
  DEFINES += -DINVISCID
endif

CPPFLAGS = -I.

# ----------- Compiler options for specific machines -------------------------

ifeq ($(ARCH),IRIX)		# -- SGI MIPS R4xxx IRIX5x
CC	 = cc
LD	 = cc
CFLAGS   = -mips2 -fullwarn $(DEFINES) -O2
LDFLAGS  = $(CFLAGS) -lmalloc -lm
endif

ifeq ($(ARCH),IRIX64)		# -- SGI MIPS R8xxx IRIX6x
CC       = cc
LD       = cc
CFLAGS   = -mips4 -O3 $(DEFINES)
LDFLAGS  = $(CFLAGS) -lmalloc -lm
endif

ifeq ($(ARCH),OSF1)		# -- DEC OSF1.
CC       = cc
LD	 = cc
CFLAGS   = -O4 -migrate -tune ev5 $(DEFINES)
LDFLAGS  = -O4 -migrate -tune ev5 $(CFLAGS) -lm
endif

ifeq ($(ARCH),Linux)		# -- GNU/Linux system.
CC	 = gcc
LD	 = gcc
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3
endif
CFLAGS   = $(DEFINES) $(OPT)
LDFLAGS  = $(CFLAGS) -lm
ifdef PROF
 CFLAGS += -p
endif
endif

# ----------------------------------------------------------------------------
#

HEADER = iso.h

all: $(PROG) IC1 IC2 IC3

# ----------------------------------------------------------------------------
# Make iso application.

ifdef TG
  ISOOBJ = allocate.o io.o nonlinear.o integrate.o energy.o FFT.o	\
	   misc.o operator.o taylor.o drive.o
else
  ISOOBJ = allocate.o io.o nonlinear.o integrate.o energy.o FFT.o	\
	   misc.o operator.o drive.o
endif

$(PROG): $(ISOOBJ) $(HEADER)
	$(LD) $(ISOOBJ) -o $(PROG) $(LDFLAGS)

$(ISOOBJ):	$(HEADER)


# ----------------------------------------------------------------------------
# Make a shaped spectrum initial conditions file.
IC1OBJ = allocate.o io.o energy.o random.o filter.o FFT.o misc.o	\
	 operator.o makeIC1.o

IC1: $(IC1OBJ) $(HEADER)
	$(LD) $(IC1OBJ) -o makeIC1 $(LDFLAGS)

# ----------------------------------------------------------------------------
# Make ICs for 3D Taylor--Green vortex.
IC2OBJ = allocate.o io.o energy.o FFT.o taylor.o misc.o operator.o makeIC2.o

IC2: $(IC2OBJ) $(IC2NNL) $(HEADER)
	$(LD) $(IC2OBJ) -o makeIC2 $(LDFLAGS)

# ----------------------------------------------------------------------------
# Make ICs for 2D Taylor vortex.
IC3OBJ = allocate.o io.o energy.o FFT.o taylor.o misc.o	\
	 operator.o makeIC3.o
IC3CHK = allocate.o io.o energy.o FFT.o taylor.o misc.o \
         operator.o taylor2D_chk.o

IC3: $(IC3OBJ) $(IC3CHK) $(HEADER)
	$(LD) $(IC3OBJ) -o makeIC3 $(LDFLAGS)
	$(LD) $(IC3CHK) -o taylor2D_chk $(LDFLAGS)

# ----------------------------------------------------------------------------
# Convert field file to tecplot (or semtex) format, with vorticity too.
TECOBJ = allocate.o io.o energy.o operator.o misc.o FFT.o iso2tec.o

iso2tec: $(TECOBJ) $(HEADER)
	$(LD) $(CFLAGS) -o iso2tec $(TECOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Make various test codes.
TESTOBJ = allocate.o energy.o operator.o FFT.o io.o	\
	  misc.o nonlinear.o random.o
TESTS   = con_chk.o div_chk.o DFT_chk.o nl_chk.o prj_chk.o

tests: $(HEADER) $(TESTOBJ) $(TESTS)
	$(LD) $(CFLAGS) -o con_chk con_chk.o  $(TESTOBJ) $(LDFLAGS)
	$(LD) $(CFLAGS) -o div_chk div_chk.o  $(TESTOBJ) $(LDFLAGS)
	$(LD) $(CFLAGS) -o dft_chk DFT_chk.o  $(TESTOBJ) $(LDFLAGS)
	$(LD) $(CFLAGS) -o nl_chk  nl_chk.o   $(TESTOBJ) $(LDFLAGS)
	$(LD) $(CFLAGS) -o prj_chk prj_chk.o  $(TESTOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Calculate the 1-D spectrum and statistics of a field dump.

SPECOBJ = allocate.o io.o energy.o FFT.o misc.o operator.o spectrum.o

spectrum: $(HEADER) $(SPECOBJ)
	$(LD) $(SPECOBJ) -o spectrum $(LDFLAGS)

# ----------------------------------------------------------------------------
# Make a source distrbution.

dist:
	tar cvf iso.tar RCS *.c *.h Makefile
	gzip iso.tar

# ----------------------------------------------------------------------------
# Clean up.
clean:
	rm -f *.o *.u *~

# ----------------------------------------------------------------------------
# Clean up.
veryclean:
	rm -f *.o *.u *~ iso makeIC1 makeIC2 makeIC3 taylor2D_chk

# ----------------------------------------------------------------------------
# S3 LES

LESOBJ = allocate.o io.o nonlin-s3.o filter.o integrate.o energy.o FFT.o \
	 misc.o operator.o drive.o

les-s3: $(LESOBJ) $(HEADER)
	$(LD) $(LESOBJ) -o les-s3 $(LDFLAGS)
