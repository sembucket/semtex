#
# mason
#

include $(PRISM)/system/defaults

PROG = mason
OBJS = mason.o load.o vertex.o edges.o nodes.o patch.o output.o \
	bandwidth.o greedy.o rcm.o sparse.o
LIBS = -lm

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROG) $(OBJS) $(LIBS)

install: mason
	-cp $(PROG) $(PRISM)/bin/$(PROG)

uninstall:
	-rm -f $(PRISM)/bin/$(PROG)

$(OBJS): mason.h

tidy:
	-rm -f *~ \#*\#
clean: tidy
	-rm -f *.o
empty: clean
	-rm -f $(PROG)

# ------------------------------------------------------------------

$(OBJS): $(PRISM)/prism/config.h
