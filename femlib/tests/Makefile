##############################################################################
# Compile test codes for femlib library functions.
#
# $Id$
##############################################################################

SEMHOME  = ../..

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(SEMHOME)/src/Makefile

# ----------------------------------------------------------------------------
# Pattern rules are enough for all codes here:
#
% : %.c
	$(CC) $(CFLAGS) -o $@ $@.c -I$(SEMHOME)/include -I. \
	-L$(SEMHOME)/lib -lfem -lvec -lm

% : %.C
	$(CXX) -o $@ $@.C $(CXXFLAGS) -I$(SEMHOME)/include -I. \
	-L$(SEMHOME)/lib/$(ARCH) -lfem -lvec $(F77LIBS) -lm
