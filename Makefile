##############################################################################
# (GNU) Makefile for spectral element solvers.
#
# $Id$
##############################################################################

tar:
	tar cvf semtex.tar *
	gzip semtex.tar

# ----------------------------------------------------------------------------
dist:
	tar cvf semtex.tar			\
	README Makefile include src lib mesh sm doc	\
	dns elliptic utility
	gzip semtex.tar

# ----------------------------------------------------------------------------
srcdist:
	tar cvf srcdist.tar				\
	README Makefile include src alplib femlib mesh	\
	dns elliptic utility
	gzip srcdist.tar

# ----------------------------------------------------------------------------
# Run this to compile libraries required by program.

libs:
	cd alplib;		\
	gmake install;		\
	gmake clean;		\
	gmake;			\
	gmake install

	cd femlib;		\
	gmake install;		\
	gmake clean;		\
	gmake;			\
	gmake install

# ----------------------------------------------------------------------------
# Run this to compile all executables.

all:
	cd src;      gmake install
	cd utility;  gmake clean; gmake all
	cd elliptic; gmake clean; gmake
	cd dns;      gmake clean; gmake

# ----------------------------------------------------------------------------
# Clean up.

clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
