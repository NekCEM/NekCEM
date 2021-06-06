prefix ?= /usr/local
bindir ?= /usr/local/bin
curdir = $(shell pwd)


build:
	@echo "Nothing to build; see README.md for details."

install:
	@mkdir -p $(DESTDIR)$(prefix)/NekCEM
	@cp -r src $(DESTDIR)$(prefix)/NekCEM
	@cp -r bin $(DESTDIR)$(prefix)/NekCEM
	@cp -r 3rd_party $(DESTDIR)$(prefix)/NekCEM
	@mkdir -p $(DESTDIR)$(bindir)
	@ln -s $(DESTDIR)$(prefix)/NekCEM/bin/configurenek $(DESTDIR)$(bindir)/configurenek
	@ln -s $(DESTDIR)$(prefix)/NekCEM/bin/nek $(DESTDIR)$(bindir)/nek

uninstall:
	@unlink $(DESTDIR)$(bindir)/configurenek
	@unlink $(DESTDIR)$(bindir)/nek
	@rm -rf $(DESTDIR)$(prefix)/NekCEM

install_inplace:
	@mkdir -p $(DESTDIR)$(bindir)
	@ln -s $(curdir)/bin/configurenek $(DESTDIR)$(bindir)/configurenek
	@ln -s $(curdir)/bin/nek $(DESTDIR)$(bindir)/nek

uninstall_inplace:
	@unlink $(DESTDIR)$(bindir)/configurenek
	@unlink $(DESTDIR)$(bindir)/nek
