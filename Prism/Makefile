#
# Prism Makefile
#
# This is the top-level Makefile, intended to be called by other applications 
# that link with various bits of Prism source code.  The variable DIM must be
# supplied as 2 or 3, and is assigned the value of 2 by default.  There are 
# four other variables that define the link between the body of Prism and 
# other applications:
#
#      APPL       application name (prism by default)
#      USER       application-specific binding (null by default)
#      SPECIAL	  application-specific source code (null by default)
#      FLAGS      application-specific compiler flags (null by default)
#
# If these variables are not overridden on the command line, then you will get
# a plain executable version of prism2d or prism3d. 
#
# There are two ways that Prism can be extended.  One is to simply provide
# a routine for performing additional analysis of the flow.  For example, to
# compile in a file "foo.c" with an internal binding to the function 
# analyze_foo(Domain *) you can simply invoke this Makefile with the args:
#
#    make USER=analyze_foo SPECIAL=foo
#
# The variable USER defines a function to be called every time step with a 
# pointer to a Domain, and SPECIAL defines a list of additional source
# files to compile and link with the main code.
#
# To create a new application you need to supply a source file with its own
# main() that coordinates execution of the time integration.  A trivial exten-
# sion of the above example with a front-end supplied in "myapp.c" would be
#
#    make USER=analyze_foo SPECIAL=foo APPL=myapp
#
# The executable is then called myapp2d or myapp3d, depdending on DIM.
#
# Finally, the variable FLAGS is used to define additional compile-time flags
# for any code supplied by the application developer.  In the above example, 
# a flag to switch on the symbols "RED" and "BLUE" could be added as
#
#    make USER=analyze_foo SPECIAL=foo APPL=myapp FLAGS="-DRED -DBLUE"
#
# Usually these variables are defined in a local Makefile that simply calls
# this one to handle compilation.  A template is provided in "template.mk".
#
# For further details see the Programmer's Guide.
# -----------------------------------------------------------------------------

include $(PRISM)/system/defaults

# These variables define the application name, binding, source code and flags

APPL    = prism
USER    =
SPECIAL =
FLAGS   =

# Dimension of the problem...must be 2 or 3

ifndef DIM
DIM = 2
endif

# Define the program to compile, source code, libraries, etc.

PROG  = $(APPL)$(DIM)d
BODY  = drive      domain    prepost   analyze    pressure  timestep \
	rotational skewsymm  stokes    fourier    flowrate  measure  \
	stats
OBJS  = $(foreach module, $(APPL) $(BODY) $(SPECIAL), $(module).o)
VPATH = $(PRISM)/prism
LIBS  = -L$(PRISM)/lib \
	$(MORTAR) $(SPECLIB) $(VECLIB) $(LAPACK) $(BLAS) $(COMM) $(SYSTEM)

# Set default and optional compilation flags

CPPFLAGS = -I$(PRISM) -I$(PRISM)/include $(OPTIONS) -DDIM=$(DIM) 

ifdef FLAGS
CPPFLAGS := $(CPPFLAGS) $(FLAGS)
endif

ifdef USER
CPPFLAGS := $(CPPFLAGS) -DUSER=$(USER)
endif

# ------------------------------------------------------------------
#               DO NOT MAKE CHANGES BELOW THIS LINE
# ------------------------------------------------------------------

$(PROG): $(OBJS)
	$(LD) $(LDFLAGS) -o $(PROG) $(OBJS) $(LIBS)

$(OBJS): $(PRISM)/prism/config.h

tidy:
	/bin/rm -f *~ \#*\#
clean: tidy
	/bin/rm -f *.o
empty: clean
	/bin/rm -f $(PROG)
