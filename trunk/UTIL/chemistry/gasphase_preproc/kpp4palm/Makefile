# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh

# ----------------------------------------------

srcdir     = ./src
bindir     = ./bin
includedir = ./include
PROG       = kpp4palm.exe

# ----------------------------------------------

all: $(PROG)

$(PROG):
	cd $(srcdir) ; \
	$(MAKE)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*~
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/$(PROG)

distclean: clean
	rm -f $(bindir)/kp4.exe

# ----------------------------------------------
