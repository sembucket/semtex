##############################################################################
# Makefile for spectral element solver utilities.
#
# $Id$
##############################################################################

SEM   = ..
VPATH = $(SEM)/src

essential: sem compare enumerate

all: essential plain nosrc meshpr calc addfield sem2nek \
     convert project interp probe rectmesh transform

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

# ----------------------------------------------------------------------------
# Build compare, exact solution comparison/restart utility.
#
COMOBJ = compare.o bsys.o BCmgr.o boundary.o condition.o feml.o mesh.o	\
	 element.o auxfield.o geometry.o pressure.o integration.o misc.o

compare: $(COMOBJ)
	$(LD) -o $@ $(COMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build enumerate, global mesh numbering utility.
#
enumerate: feml.o mesh.o enumerate.o
	$(LD) -o $@ $@.o feml.o mesh.o $(LDFLAGS)

# ----------------------------------------------------------------------------
# Plain (C) programmes.
#
plain: slit chop
slit chop:\
	slit.c chop.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.c -lm

# ----------------------------------------------------------------------------
# Standalone (bar semtex standard libs) C programs.
#
nosrc: avgdump moden noiz rstress xplane sem2tec convert
avgdump moden noiz rstress xplane sem2tec convert: \
	avgdump.c moden.c xplane.c noiz.c rstress.c sem2tec.c convert.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $@.c $(LDFLAGS) -lm

# ----------------------------------------------------------------------------
# Build meshpr, 2D mesh node computation utility.
#
meshpr:	feml.o mesh.o meshpr.o
	$(LD) -o $@ $@.o feml.o mesh.o $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build calc, an interactive interface to the alplib function parser.
#
calc: calc.o
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build addfield, field postprocessing utility.
#
ADFOBJ = addfield.o auxfield.o BCmgr.o boundary.o bsys.o condition.o 	\
	 domain.o element.o feml.o field.o geometry.o integration.o	\
	 mesh.o misc.o pressure.o

addfield: $(ADFOBJ)
	$(LD) -o $@ $(ADFOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build addvort utility (just compute the vorticity field).
#
ADDOBJ = addvort.o auxfield.o BCmgr.o boundary.o bsys.o condition.o 	\
	 domain.o element.o feml.o field.o geometry.o integration.o	\
	 mesh.o misc.o pressure.o

addvort: $(ADDOBJ)
	$(LD) -o $@ $(ADDOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build addvortr utility, line addvort but creates another field x^2+y^2.
#
ADROBJ = addvortr.o auxfield.o BCmgr.o boundary.o domain.o	\
	 element.o feml.o field.o geometry.o  integration.o mesh.o	\
	 misc.o	pressure.o \

addvortr: $(ADROBJ)
	$(LD) -o $@ $(ADROBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build semtex --> NEKTON input file converter.
#
sem2nek: sem2nek.o feml.o mesh.o
	$(LD) -o $@ $@.o feml.o mesh.o $(S2NOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build project, field dump projection to higher/lower orders.
#
project: project.C
	$(CXX) $(CPPFLAGS) -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build transform, which does forward/inverse Fourier/Legendre transforms.
#
transform: transform.C
	$(CXX) $(CPPFLAGS) -g -o $@ $@.C $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build interp, field dump interpolator.
#
INTOBJ = interp.o feml.o mesh.o element.o auxfield.o geometry.o

interp: $(INTOBJ)
	$(CXX) $(CPPFLAGS) -o $@ $(INTOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build integral, field dump scalar integral evaluator.
#
integral: integral.o  $(INTOBJ)
	$(CXX) $(CPPFLAGS) -o $@ $(INTOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build probe, field dump data extraction.
#
PRBOBJ = probe.o feml.o mesh.o element.o auxfield.o geometry.o

probe: $(PRBOBJ)
	$(CXX) $(CPPFLAGS) -o $@ $(PRBOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build rectmesh, rectangular mesh generator.
#
rectmesh: rectmesh.C
	$(CXX) $(CPPFLAGS) -o $@ $@.C -g

# ----------------------------------------------------------------------------
# Special machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.C
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.C
endif
