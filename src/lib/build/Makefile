# Makefile
#
# This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
#
# COPYRIGHT (c) 2016 The University of Chicago
# All rights reserved.
#
# The Makefile for building all versions of the Diderot runtime system
#
# src/lib/build/Makefile.  Generated from Makefile.in by configure.
#

SHELL =			/bin/bash
INSTALL =		/usr/bin/install -c


# clear out unused file suffixes for implicit rules
#
.SUFFIXES:
.SUFFIXES:	.cxx .hxx

#
# directories for the install target
#
PREFIX =		/usr/local
INSTALL_LIBDIR =	$(PREFIX)/lib

#
# directories for the local-install target
#
LIBDIR =	/home/teocollin/gitcode/diderot/lib

#
# optional logging support
#
ifeq (false,true)
LOGGING_HXX =		/home/teocollin/gitcode/diderot/src/lib/include/diderot/logging.hxx
endif

.PHONY:			build local-install install

build:			$(LOGGING_HXX)
	for dir in $(wildcard rt-*) ; do \
	  (cd $$dir; $(MAKE) build) || exit $$?; \
	done

local-install:		$(LOGGING_HXX)
	$(INSTALL) -d $(LIBDIR)
	for dir in $(wildcard rt-*); do \
	  (cd $$dir && $(MAKE) local-install) || exit $$?; \
	done

install:		$(LOGGING_HXX)
	$(INSTALL) -d $(INSTALL_LIBDIR)
	for dir in $(wildcard rt-*); do \
	  (cd $$dir && $(MAKE) install) || exit $$?; \
	done

#
# optional logging support; note that this target is redundant with the
# logging target in the top-level Makefile
#
$(LOGGING_HXX):
	(cd /home/teocollin/gitcode/diderot/src/logging/gen; log-gen.sh) || exit $$?

#################### Cleanup ####################

CLEAN_SUBDIRS =		$(wildcard rt-*)
DISTCLEAN_FILES =	$(wildcard rt-*) mk/build.gmk
DEVCLEAN_FILES =

include /home/teocollin/gitcode/diderot/mk/clean-rules.gmk
