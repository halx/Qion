# $Id: Makefile,v 1.5 2005/06/07 10:12:14 hal Exp $
#
# top-level Makefile for the QM/MM-MD program
#



SHELL	= /bin/sh
RM	= rm -f

# program name
PROGRAM	= Qion

# directories
SRCDIR = source
DOCDIR = documents
TESTDIR = testdir

SUBDIR = $(SRCDIR) $(DOCDIR)

# common flags for GNU compilers
GNUFLAGS = -g3 -Wall -W

# GNU Fortran 77
FC	= gfortran
FFLAGS	= $(GNUFLAGS) #-fpedantic

# PGI compiler
#FC	= pgf77
#FFLAGS	= -fast -Mdclchk -Mvect=assoc,cachesize:524288 #-Mconcur

# Intel Compiler (for highest performance)
#FC      = ifort
#FFLAGS  = -ipo -O3 -unroll -FI -W0 -i_dynamic
#LDFLAGS	= -Vaxlib

# Compaq Fortran compiler
#FC	= fort
#FFLAGS	= -O5 -fast

# AIX xlf
#FC	= f77
#FFLAGS	= -O3 -qextname

# C compiler
CC	= gcc
CFLAGS	= $(GNUFLAGS)

#CC	= pgcc
#CFLAGS	= -O2

#LDFLAGS	+= -s



all: prog

.PHONY: prog docs clean testclean distclean


prog:
	cd ./$(SRCDIR); $(MAKE) FC='$(FC)' FFLAGS='$(FFLAGS)' CC='$(CC)' \
		CFLAGS='$(CFLAGS)' LDFLAGS='$(LDFLAGS)' CPPFLAGS='$(CPPFLAGS)' \
		PROGRAM='../$(PROGRAM)'

docs:
	cd ./$(DOCDIR); $(MAKE)


clean:
	for d in $(SUBDIR); do (cd ./$$d && $(MAKE) $@); done
	$(RM) core

testclean:
	cd $(TESTDIR); $(RM) qion.rst qion.out qion.info qion.en qion.traj\
	qion.vel qion_umb.hst tm.out alpha basis beta mos control coord \
	energy gradient	ddens dens* fock natural statistics test.* gauss.in \
	gauss.out qmmm_* *.prev *~ core

distclean: testclean
	for d in $(SUBDIR); do (cd ./$$d && $(MAKE) $@); done
	$(RM) $(PROGRAM)
