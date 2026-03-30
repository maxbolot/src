FC       = ifx
SHELL    = /bin/bash
FCFLAGS  = -O2 -fPIC
MODULE_INIT = source /usr/share/Modules/init/bash >/dev/null 2>&1 || source /etc/profile.d/modules.sh >/dev/null 2>&1
MODULE_LOAD = module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

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
	$(MODULE_INIT) && \
	$(MODULE_LOAD) && \
	$(FC) $(FCFLAGS) -c $< -o $@

install: $(LIB)
	install -d $(LIBDIR) $(INCDIR)
	install -m 644 $(LIB) $(LIBDIR)/
	install -m 644 $(MODS) $(INCDIR)/

uninstall:
	rm -f $(LIBDIR)/$(LIB) $(INCDIR)/$(MODS)

clean:
	rm -f $(OBJS) $(MODS) $(LIB)
