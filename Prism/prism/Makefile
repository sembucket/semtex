#
# Generic Recursive Makefile
#
# This makefile has six targets, listed below, to perform the 
# basic operation required for managing a collection of libraries.
# Each target loops over all subdirectories with a Makefile and
# execute Make on them with the same target.
# 
# $Id$

MODULE = prism
INCDIR = $(PRISM)/include
LIBDIR = 

all:
	for dir in `ls`; do \
		if [ -f $$dir/Makefile ]; \
		then ( cd $$dir; $(MAKE)); \
		fi; \
	done

install:
	for dir in `ls`; do \
		if [ -f $$dir/Makefile ]; \
		then ( cd $$dir; $(MAKE) install); \
		fi; \
	done 
	if [ ! -d $(INCDIR)/$(MODULE) ]; \
	then mkdir $(INCDIR)/$(MODULE); \
	else $(RM) $(INCDIR)/$(MODULE)/* ; \
	fi
	-cp *.h $(INCDIR)/$(MODULE)/. 

uninstall: empty
	for dir in `ls`; do \
		if [ -f $$dir/Makefile ]; \
		then ( cd $$dir; $(MAKE) uninstall); \
		fi; \
	done
	$(RM) -r $(INCDIR)/$(MODULE)

tidy:
	for dir in `ls`; do \
		if [ -f $$dir/Makefile ]; \
		then ( cd $$dir; $(MAKE) tidy); \
		fi; \
	done


clean: tidy
	for dir in `ls`; do \
		if [ -f $$dir/Makefile ]; \
		then ( cd $$dir; $(MAKE) clean); \
		fi; \
	done

empty: clean
	for dir in `ls`; do \
		if [ -f $$dir/Makefile ]; \
		then ( cd $$dir; $(MAKE) empty); \
		fi; \
	done


