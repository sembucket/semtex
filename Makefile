##############################################################################
# (GNU) Makefile for spectral element solvers.
#
# $Id$
##############################################################################

# -- MAKE supplies the path to GNU make.

MAKE = gmake

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
	tar cvf srcdist.tar				  \
	README Makefile include src alplib femlib mesh sm \
	dns elliptic utility
	gzip srcdist.tar

# ----------------------------------------------------------------------------
# Run this to compile libraries required by program.

libs:
	cd alplib;		\
	$(MAKE) -s install;	\

	cd femlib;		\
	$(MAKE) -s install;	\

	cd alplib;		\
	$(MAKE) -s clean;	\
	$(MAKE) -s;		\
	$(MAKE) -s install

	cd femlib;		\
	$(MAKE) -s clean;	\
	$(MAKE) -s;		\
	$(MAKE) -s install

# ----------------------------------------------------------------------------
# Make version of femlib with MPI.

parlib:

	cd femlib;		\
	$(MAKE) -s install;	\

	cd femlib;		\
	$(MAKE) -s clean;	\
	$(MAKE) -s MPI=1;	\
	$(MAKE) -s install MPI=1


# ----------------------------------------------------------------------------
# Run this to compile all executables.

all:
	cd src;      $(MAKE) install
	cd utility;  $(MAKE) clean; $(MAKE) all
	cd elliptic; $(MAKE) clean; $(MAKE)
	cd dns;      $(MAKE) clean; $(MAKE)

# ----------------------------------------------------------------------------
# Run a test of the (serial) DNS solver.  This could take a few minutes.

test:  libs
	cd utility; $(MAKE) -s clean; $(MAKE) -s enumerate; $(MAKE) -s compare
	cd dns; $(MAKE) -s clean; $(MAKE) -s ;		\
	rm -f compare;   ln -s ../utility/compare   . ;	\
	rm -f enumerate; ln -s ../utility/enumerate .
	@echo -- No output from testregress indicates success. --
	cd dns; testregress dns

# ----------------------------------------------------------------------------
# Run test of parallel version of DNS solver: do "make test" first.
# Also, you may need to edit the file dns/testregress_mp to get MPI to run.

partest: parlib
	cd dns; $(MAKE) -s clean; $(MAKE) -s ALIAS=1 MPI=1;
	@echo -- No output from testregress_mp indicates success. --
	cd dns; testregress_mp dns_mp

# ----------------------------------------------------------------------------
# Clean up.

clean:
	rm -f *.o *~
	rm -rf ILDUMPS
	rm -rf ii_files
