FC       = ifx
FCFLAGS  = -O2 -fPIC

LIB      = libhist.a
OBJS     = hist.o
MODS     = hist.mod

DESTDIR  = $(HOME)/local
LIBDIR   = $(DESTDIR)/lib
INCDIR   = $(DESTDIR)/include

.PHONY: all install uninstall clean

all: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

hist.o: hist.f90
	$(FC) $(FCFLAGS) -c $< -o $@

install: $(LIB)
	install -d $(LIBDIR) $(INCDIR)
	install -m 644 $(LIB) $(LIBDIR)/
	install -m 644 $(MODS) $(INCDIR)/

uninstall:
	rm -f $(LIBDIR)/$(LIB) $(INCDIR)/$(MODS)

clean:
	rm -f $(OBJS) $(MODS) $(LIB)
