# ----------------------------------------------------------------------------
# Define architectures.
#

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

RECIPES = -L$(HOME)/lib/$(ARCH) -L/usr/local/lib -lrecipesDP
SEMHOME = ../semtex
SEMHDR  = -I$(SEMHOME)/include
VECLIB  = -L$(SEMHOME)/lib/$(ARCH)/ -lalp

###############################################################################

asp: FFTutil.o asp.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o -lm

acor: FFTutil.o acor.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o -lm

correlate: FFTutil.o correlate.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o $(RECIPES) -lm

chop: chop.c
	$(CC) $(CFLAGS) -o chop chop.c

slit: FFTutil.o slit.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o -lm

pdf: FFTutil.o pdf.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o -lm

memspec: FFTutil.o memspec.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o -lm

dft_chk: FFTutil.o dft_chk.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o -lm

pcor: pcor.c
	$(CC) $(CFLAGS) -o $@ $@.c

upxf: upxf.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

peaks: peaks.C
	$(CXX) -o $@ $@.C

cycint: cycint.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

cycle: cycle.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

ratfit: ratfit.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

smooth: smooth.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

embed: embed.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

odeint: odeint.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

zeros: zeros.C
	$(CXX) -o $@ $@.C $(CXXFLAGS)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(SEMHOME)/src/Makefile

