FC       = ifx
FCFLAGS  = -O2 -fPIC
NCFLAGS  = $(shell nf-config --fflags)

LIB      = libnc.a
OBJS     = nc.o
MODS     = nc.mod

DESTDIR  = $(HOME)/local
LIBDIR   = $(DESTDIR)/lib
INCDIR   = $(DESTDIR)/include

.PHONY: all install uninstall clean

all: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

nc.o: nc.f90
	$(FC) $(FCFLAGS) $(NCFLAGS) -c $< -o $@

install: $(LIB)
	install -d $(LIBDIR) $(INCDIR)
	install -m 644 $(LIB) $(LIBDIR)/
	install -m 644 $(MODS) $(INCDIR)/

uninstall:
	rm -f $(LIBDIR)/$(LIB) $(INCDIR)/$(MODS)

clean:
	rm -f $(OBJS) $(MODS) $(LIB)
