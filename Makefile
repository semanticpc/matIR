-include MakeDefns

INSTALLDIRS = $(bindir) $(includedir) $(pkgincludedir)  $(libdir) 

.PHONY: all dist clean install $(INSTALLDIRS) site-search

all: 
	$(MAKE) -C obj -f ../src/Makefile
	$(MAKE) -C utils
	

$(INSTALLDIRS):
	$(INSTALL_DIR) $@

clean:
	$(MAKE) clean -C utils
	$(MAKE) clean -C obj -f ../src/Makefile
	rm -f depend/*

install: $(INSTALLDIRS)
	rm -f $(libdir)/$(INDRILIB)
	$(MAKE) install -C utils
	$(MAKE) install -C obj -f ../src/Makefile
	$(INSTALL_DATA) Makefile.app $(pkgdatadir)

test:
