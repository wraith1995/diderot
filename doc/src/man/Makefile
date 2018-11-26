# Makefile for diderot man pages
#
# This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
#
# COPYRIGHT (c) 2017 The University of Chicago
# All rights reserved.
#

SHELL =		/bin/bash


ASCIIDOCTOR =	none
ATTRIBUTES =	-a version=1.0

# set attributes for conditional parts of the documentation
#
ifeq (false,true)
ATTRIBUTES +=	-a debugger_enabled
endif
ifeq (false,true)
ATTRIBUTES +=	-a opencl_enabled
endif
ifeq (false,true)
ATTRIBUTES +=	-a cuda_enabled
endif
ATTRIBUTES +=	-a teem_dir="/home/teocollin/installs/teem/teem-build"

MANPAGES =	diderotc.1 diderot-api.5
HTML =		$(patsubst %,%.html,$(MANPAGES))

.PHONY:		all
all:		$(MANPAGES) $(HTML)

diderotc.1:	diderotc.1.txt
	$(ASCIIDOCTOR) $(ATTRIBUTES) -b manpage diderotc.1.txt

diderot-api.5:  diderot-api.5.txt
	$(ASCIIDOCTOR) $(ATTRIBUTES) -b manpage diderot-api.5.txt

diderotc.1.html:	diderotc.1.txt
	$(ASCIIDOCTOR) $(ATTRIBUTES) -b html5 diderotc.1.txt

diderot-api.5.html:  diderot-api.5.txt
	$(ASCIIDOCTOR) $(ATTRIBUTES) -b html5 diderot-api.5.txt

#################### Cleanup ####################

CLEAN_SUBDIRS =

CLEAN_FILES =		$(MANPAGES) $(HTML)

DISTCLEAN_FILES +=	Makefile

DEVCLEAN_FILES =

include /home/teocollin/gitcode/diderot/mk/clean-rules.gmk
