#
# Template Makefile for Prism-based applications
#

include $(PRISM)/system/defaults

# Define each of the following variables!  The first three must define a
# single symbol (DIM=2 APPL=myapp USER=myinterface), but the last two can be a 
# list of symbols (FLAGS=-DA -DB -DC).   Note that source files defined in 
# SPECIAL should be given w/out extensions (SPECIAL=file1 file2 file3).

DIM     =
APPL    = 
USER    = 
SPECIAL = 
FLAGS   =

# Now invoke the Prism makefile to handle compilation

$(APPL):
	$(MAKE) -f $(PRISM)/Makefile DIM=$(DIM) \
	APPL=$(APPL) USER=$(USER) SPECIAL="$(SPECIAL)" FLAGS="$(FLAGS)"

tidy:
	-rm -f *~ \#*\#
clean: tidy
	-rm -f *.o core
empty: clean
	-rm -f $(APPL)*d

