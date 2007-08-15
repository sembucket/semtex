##############################################################################
# Makefile for spectral element solver utilities.
#
# $Id$
##############################################################################

SEM   = ..
VPATH = $(SEM)/src

essential: sem compare enumerate

all: essential plain preplot nosrc meshpr calc addfield sem2nek \
     convert project interp probe rectmesh transform wallmesh   \
     phase lowpass eneq integral traction

distclean:
	make clean
	rm -f addfield avgdump calc chop compare convert eneq enumerate \
        integral interp lowpass meshpr moden noiz phase preplot probe   \
        probeline probeplane project rectmesh repeatxy repeatz rstress  \
        sem2nek sem2tec slit transform wallmesh xplane

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile
SEMHDR = sem.h

# ----------------------------------------------------------------------------
# Build compare, exact solution comparison/restart utility.
#
COMOBJ = compare.o bsys.o bcmgr.o boundary.o condition.o edge.o		\
	 family.o feml.o mesh.o element.o auxfield.o geometry.o		\
         pressure.o integration.o misc.o
$(COMOBJ): $(SEMHDR)

compare: $(COMOBJ)
	$(LD) -o $@ $(COMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build enumerate, global mesh numbering utility.
#
enumerate: feml.o mesh.o enumerate.o $(SEMHDR)
	$(LD) -o $@ $@.o feml.o mesh.o $(LDFLAGS)

# ----------------------------------------------------------------------------
# Plain (C) programmes.
#

# -- Default build for a C source file.

.c: 
	$(CC) -o $@ $@.c $(CFLAGS) $(CPPFLAGS) $(CLDFLAGS) -lm

plain: slit chop
slit chop:\
	slit.c chop.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.c -lm

preplot: preplot.c
	$(CC) -DIRIS $(CFLAGS) $(CPPFLAGS) -o $@ $@.c -lm

# ----------------------------------------------------------------------------
# Standalone (bar semtex standard libs) C programs.
#

nosrc: avgdump moden noiz rstress xplane sem2tec convert repeatz repeatxy
avgdump moden noiz rstress xplane sem2tec convert repeatz repeatxy : \
	avgdump.c moden.c xplane.c noiz.c rstress.c sem2tec.c \
        convert.c repeatz.c repeatxy.c
	$(CC) -o $@ $@.c $(CFLAGS) $(CPPFLAGS) $(CLDFLAGS) -lm

# ----------------------------------------------------------------------------
# Build meshpr, 2D mesh node computation utility.
#
meshpr:	feml.o mesh.o meshpr.o $(SEMHDR)
	$(LD) -o $@ $@.o feml.o mesh.o $(LDFLAGS)


# ----------------------------------------------------------------------------
# Build wallmesh, 2D mesh filter utility for use with meshpr.
#
WALOBJ = wallmesh.o bsys.o bcmgr.o condition.o pressure.o edge.o	\
	 family.o feml.o mesh.o element.o geometry.o integration.o
         
$(WALOBJ): $(SEMHDR)

wallmesh: $(WALOBJ)
	$(LD) -o $@ $(WALOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build traction, utility to compute normal & tangential tractions on walls
#
TRCOBJ = traction.o feml.o auxfield.o bcmgr.o element.o geometry.o	\
         domain.o condition.o family.o mesh.o bsys.o field.o edge.o	\
         boundary.o pressure.o misc.o integration.o
         
$(TRCOBJ): $(SEMHDR)

traction: $(TRCOBJ)
	$(LD) -o $@ $(TRCOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build calc, an interactive interface to the femlib function parser.
#
calc: calc.o
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build addfield, field postprocessing utility.
#
ADFOBJ = addfield.o auxfield.o bcmgr.o boundary.o bsys.o condition.o 	\
	 domain.o edge.o element.o family.o feml.o field.o geometry.o	\
	 integration.o mesh.o misc.o pressure.o tensorcalcs.o
$(ADFOBJ): $(SEMHDR) tensorcalcs.h

addfield: $(ADFOBJ)
	$(LD) -o $@ $(ADFOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build variant of addfield that multiplies enstrophy by sign of another var.
#
SIGOBJ = signens.o auxfield.o bcmgr.o boundary.o bsys.o condition.o 	\
	 domain.o edge.o element.o family.o feml.o field.o geometry.o	\
	 integration.o mesh.o misc.o pressure.o
$(SIGOBJ): $(SEMHDR)

signens: $(SIGOBJ)
	$(LD) -o $@ $(SIGOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build semtex --> NEKTON input file converter.
#
sem2nek: sem2nek.o feml.o mesh.o $(SEMHDR)
	$(LD) -o $@ $@.o feml.o mesh.o $(S2NOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build project, field dump projection to higher/lower orders.
#
project: project.C
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build transform, which does forward/inverse Fourier/Legendre transforms.
#
transform: transform.o data2df.o data2df.h $(SEMHDR)
	$(CXX) $(CPPFLAGS) -o $@ $@.o data2df.o $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build phase, which performs operations on 3D data in phase/Fourier space.
#
phase: phase.o data2df.o data2df.h $(SEMHDR)
	$(CXX) $(CPPFLAGS) -o $@ $@.o data2df.o $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build lowpass, which filters data.
#
lowpass: lowpass.o data2df.o data2df.h $(SEMHDR)
	$(CXX) $(CPPFLAGS) -o $@ $@.o data2df.o $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build interp, field dump interpolator.
#
INTOBJ = feml.o mesh.o element.o family.o auxfield.o geometry.o
$(INTOBJ): $(SEMHDR)

interp: interp.o $(INTOBJ)
	$(CXX) $(CPPFLAGS) -o $@ $@.o $(INTOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build integral, field dump scalar integral evaluator.
#
integral: integral.o  $(INTOBJ)
	$(CXX) $(CPPFLAGS) -o $@ $@.o $(INTOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build probe, field dump data extraction.
#
PRBOBJ = probe.o feml.o mesh.o element.o family.o auxfield.o geometry.o
$(PRBOBJ): $(SEMHDR)

probe: $(PRBOBJ)
	$(CXX) $(CPPFLAGS) -o $@ $(PRBOBJ) $(LDFLAGS)
	ln -s -f probe probeline
	ln -s -f probe probeplane

# ----------------------------------------------------------------------------
# Build eneq, which computes terms in the fluctuating flow energy equations.
#
ENOBJ = eneq.o feml.o mesh.o element.o family.o auxfield.o geometry.o misc.o
$(ENOBJ): $(SEMHDR)

eneq: $(ENOBJ)
	$(CXX) $(CPPFLAGS) -o $@ $(ENOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build rectmesh, rectangular mesh generator.
#
rectmesh: rectmesh.C
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Special machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif
