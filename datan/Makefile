# ----------------------------------------------------------------------------
# Define architectures.
#
ARCH = $(shell uname -s)
MACH = $(shell uname -m)
ifeq ($(ARCH),Linux)
  ifeq ($(MACH),alpha)
    ARCH = Linux-alpha
  endif
endif

RECIPES = -L$(HOME)/lib/$(ARCH) -L/usr/local/lib -lrecipesDP
SEMHOME = ../semtex
SEMHDR  = -I$(SEMHOME)/include
VECLIB  = -L$(SEMHOME)/lib/$(ARCH)/ -lalp

# ----------------------------------------------------------------------------
# Compiler options for specific machines.
#

ifeq ($(ARCH),SunOS)              	# -- Sun.
CXX    = CC
LD     = CC
CC     = cc
ifdef DEBUG
  OPT    = -g -xarch=v9
else
  OPT    = -fast -xarch=v9
endif
CFLAGS   = $(DEFINES) $(OPT)
CXXFLAGS = $(DEFINES) $(OPT) -library=iostream
LDFLAGS  = $(LIBDIR)  $(OPT) -library=iostream
CLDFLAGS = $(LIBDIR)  $(OPT)
endif

ifeq ($(ARCH),Linux-alpha)              # -- Linux alpha w Compaq C & f90.
CXX    = g++
LD     = g++
CC     = ccc
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3
endif
#CFLAGS   = $(DEFINES) $(OPT)
CFLAGS   = $(DEFINES) -arch host -tune host -fast -O4
CXXFLAGS = $(DEFINES) $(OPT)
LDFLAGS  = $(LIBDIR)
CLDFLAGS = $(LDFLAGS)
F77LIBS  = -lfor -lots -lg2c
endif

ifeq ($(ARCH),Linux)               # -- Linux Intel w GNU compilation.
CXX      = g++
LD       = g++
CC       = gcc
ifdef DEBUG
  OPT    = -gstabs+
else
  OPT    = -O3 -w
endif
CFLAGS   = $(DEFINES) $(OPT)
CXXFLAGS = $(DEFINES) $(OPT)
CLDFLAGS = $(LDFLAGS)
F77LIBS  = -lg2c
endif

ifeq ($(ARCH),IRIX)			# -- SGI 32-bit systems.
CXX	 = CC
LD	 = C
RANLIB   = 
CFLAGS   = $(DEFINES) -fullwarn -O2
CXXFLAGS = $(DEFINES) -fullwarn -O2
LDFLAGS  = $(LIBDIR) -lfem -lalp $(F77LIBS) -lmalloc -lm
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),IRIX64)			# -- SGI 64-bit systems.
CXX	 = CC
LD	 = CC
CC       = cc
CFLAGS   = $(DEFINES) -n32 -mips4  $(OPT)
CXXFLAGS = $(DEFINES) -n32 -mips4  $(OPT)
F77LIBS  = -lftn
ifdef DEBUG
  OPT    = -g
else
  OPT    =  -O2
endif
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),OSF1)			# -- DEC OSF1 Alpha.
CXX	 = cxx
LD	 = cxx
CC       = cc
CFLAGS   = $(DEFINES) -fast -inline speed -fp_reorder \
           -speculate all -tune host -O4
CXXFLAGS = $(DEFINES) -tune host -O2
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),UNIX_System_V)    	# -- Fujitsu VPP300.
CC       = vcc
CXX      = /home/home01/565/hmb565/C++/bin/CC
LD       = $(CXX)
OPT      = -Wv,-m3,-Om,-te -Kpopt
CFLAGS   = $(DEFINES) -Wv,-m3,-Om -Kpopt
CXXFLAGS = $(DEFINES)  $(OPT)
LDFLAGS  = $(OPT) $(CPPFLAGS) $(LIBDIR) $(LIBS)
CLDFLAGS = $(LDFLAGS)
endif

ifeq ($(ARCH),SUPER-UX)         	# -- NEC SX-4.
DEFINES += -D_BUILTIN_ -D_SX
CXX       = c++
LD        = c++
CFLAGS    = $(DEFINES) -pvctl,fullmsg,loopcnt=20000 -Onooverlap \
           -hfloat0 -hnostkchk
CXXFLAGS  = -V $(DEFINES) -pvctl,loopcnt=20000,fullmsg -acct $(OPT)
CPPFLAGS += -I/usr/include/CC
OPT       = -O nooverlap -Nstkchk
NOVEC     = -c -Nstkchk -acct -dir novec
LDFLAGS   =  -verbose -V -p -f90libew $(CXXFLAGS) $(CPPFLAGS) $(LIBDIR) $(LIBS)
CLDFLAGS  = $(LDFLAGS)
endif

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
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o-lm

memspec: FFTutil.o memspec.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o -lm

dft_chk: FFTutil.o dft_chk.o $(INCLUDE)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.o FFTutil.o -lm

pcor: pcor.c
	$(CC) $(CFLAGS) -o $@ $@.c

upxf: upxf.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(RECIPES) $(F77LIBS)

peaks: peaks.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(RECIPES) $(F77LIBS)

cycint: cycint.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

cycle: cycle.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

smooth: smooth.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) $(SEMHDR) $(VECLIB) $(RECIPES) $(F77LIBS)

clean:
	rm -rf *.o *~ core
