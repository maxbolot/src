FC       = ifx
FCFLAGS  = -O2 -fPIC

LIB      = libinterp.a
OBJS     = interp.o
MODS     = interp.mod

DESTDIR  = $(HOME)/local
LIBDIR   = $(DESTDIR)/lib
INCDIR   = $(DESTDIR)/include

.PHONY: all install uninstall clean

all: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

interp.o: interp.f90
	$(FC) $(FCFLAGS) -c $< -o $@

install: $(LIB)
	install -d $(LIBDIR) $(INCDIR)
	install -m 644 $(LIB) $(LIBDIR)/
	install -m 644 $(MODS) $(INCDIR)/

uninstall:
	rm -f $(LIBDIR)/$(LIB) $(INCDIR)/$(MODS)

clean:
	rm -f $(OBJS) $(MODS) $(LIB)
