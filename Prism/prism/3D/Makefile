#
# Makefile for the 3D version of Prism
#

include $(PRISM)/system/defaults

# Define each of the following variables!  The first three must define a
# single symbol (DIM=2 APPL=myapp USER=myinterface), but the last two can be a 
# list of symbols (FLAGS=-DA -DB -DC).   Note that source files defined in 
# SPECIAL should be given w/out extensions (SPECIAL=file1 file2 file3).

DIM     = 3
APPL    = prism
USER    = 
SPECIAL = 
FLAGS   =

# Now invoke the Prism makefile to handle compilation

$(APPL)$(DIM)d:
	$(MAKE) -f $(PRISM)/Makefile DIM=$(DIM) \
	APPL=$(APPL) USER=$(USER) SPECIAL="$(SPECIAL)" FLAGS="$(FLAGS)"

install: $(APPL)$(DIM)d
ifdef PARALLEL
	cp $(APPL)$(DIM)d $(PRISM)/bin/$(APPL)$(DIM)d_mp
else
	cp $(APPL)$(DIM)d $(PRISM)/bin/$(APPL)$(DIM)d
endif

uninstall: 
	$(RM) $(PRISM)/bin/$(APPL)$(DIM)d


tidy:
	-rm -f *~ \#*\#
clean: tidy
	-rm -f *.o core
empty: clean
	-rm -f $(APPL)$(DIM)d

