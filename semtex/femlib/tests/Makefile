##############################################################################
# Compile test codes for femlib library functions.
#
# $Id$
##############################################################################

SEMHOME  = ../..
LIB      = libfem.a

ARCH = $(shell uname -s)
MACH = $(shell uname -m)

ifeq ($(ARCH),Linux)
  ifeq ($(MACH),alpha)
    ARCH = Linux-alpha
  endif
  ifeq ($(MACH),i586)
    ARCH = Linux-iX86
  endif
  ifeq ($(MACH),i686)
    ARCH = Linux-iX86
  endif
endif

DEFINES  =
ifdef DEBUG
  DEFINES += -DDEBUG
endif

CC	 = cc
CXX	 = CC
FC	 = f77
OPTIONS  = -g
CFLAGS	 = $(OPTIONS)
FFLAGS   = $(OPTIONS)
ARFLAGS	 = cr
CPPFLAGS = -I$(SEMHOME)/include
RANLIB   =

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(SEMHOME)/src/Makefile

testpoly: testpoly.c
	$(CC) $(CFLAGS) -o $@ $@.c -I$(SEMHOME)/include -I. \
	-L$(SEMHOME)/lib -lfem -lalp -lm

testFFT: testFFT.C
	$(CXX) -o $@ $@.C -I$(SEMHOME)/include -I. -g \
	$(CXXFLAGS) -L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

testFFTW: testFFTW.C
	$(CXX) -o $@ $@.C -I$(SEMHOME)/include -I. -g \
	$(CXXFLAGS) -L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

testDLT: testDLT.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

testDLT2: testDLT2.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

testDLT3: testDLT3.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

testBVD: testBVD.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

glzw: glzw.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

shapes: shapes.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

lagint: lagint.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

lagpro: lagpro.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm

dpt: dpt.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. -g \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lalp $(F77LIBS) -lm
