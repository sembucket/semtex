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

#RECIPES = -L$(HOME)/lib/$(ARCH) -L/usr/local/lib -lrecipesDP
RECIPES =  -L/sw/lib -L/sw/lib/gcc4.8/lib -lrecipesDP
SEMHOME = ../semtex
SEMHDR  = -I$(SEMHOME)/include
VECLIB  = -L$(SEMHOME)/lib/$(ARCH)/ -lvec

###############################################################################

asp: FFTutil.o asp.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o -lm

acor: FFTutil.o acor.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o $(F77LIBS) -lm

correlate: FFTutil.o correlate.o recipes.o $(INCLUDE)
	$(CC) -o $@ $@.o FFTutil.o recipes.o $(F77LIBS) -lm

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

upxf: upxf.o recipes.o
	$(CXX) -o $@ $@.o recipes.o $(VECLIB) $(F77LIBS) -lm

peaks: peaks.C
	$(CXX) -o $@ $@.C

cycint: cycint.o recipes.o
	$(CXX) -o $@ $@.o recipes.o $(VECLIB) $(F77LIBS) -lm

cycle: cycle.o recipes.o
	$(CXX) -o $@ $@.o recipes.o $(VECLIB) $(F77LIBS) -lm

ratfit: ratfit.o recipes.o
	$(CXX) -o $@ $@.o recipes.o $(VECLIB) $(F77LIBS) -lm

smooth: smooth.o recipes.o
	$(CXX) -o $@ $@.o recipes.o $(VECLIB) $(F77LIBS) -lm

embed: embed.o recipes.o
	$(CXX) -o $@ $@.o recipes.o $(VECLIB) $(F77LIBS) -lm

odeint: odeint.o recipes.o
	$(CXX) -o $@ $@.o recipes.o $(VECLIB) $(F77LIBS) -lm

zeros: zeros.C
	$(CXX) -o $@ $@.C $(CXXFLAGS)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(SEMHOME)/src/Makefile
ifeq ($(ARCH),Darwin)
FFLAGS  = -m64 -O2 -ffixed-line-length-none -fno-range-check -fno-align-commons
F77LIBS += -L/sw/lib  -L/sw/lib/gcc4.8/lib
CXXFLAGS += $(SEMHDR) -Wno-deprecated
endif
