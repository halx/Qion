# $Id: Makefile,v 1.1.1.1 2003/06/17 10:15:22 hal Exp $
#
# Makefile for the documentation of the QM/MM-MD program
#



SHELL	= /bin/sh
RM	= rm -f

# name of the document
MANUAL = manual

ABBREV = journal-abbrev.bib

# programs
LATEX	= latex
BIBTEX	= bibtex
DVIPS	= dvips
PDFLATEX = pdflatex



all: dvi ps pdf

.PHONY: dvi ps pdf clean distclean


dvi: $(MANUAL).dvi

ps: $(MANUAL).ps

pdf: $(MANUAL).pdf


$(MANUAL).dvi: $(MANUAL).tex $(MANUAL).bib $(MANUAL).bbl
	$(LATEX) $(MANUAL).tex
	$(LATEX) $(MANUAL).tex

$(MANUAL).ps: $(MANUAL).dvi
	$(DVIPS) -Ppdf -o $(MANUAL).ps $(MANUAL).dvi

$(MANUAL).pdf: $(MANUAL).dvi
	$(PDFLATEX) $(MANUAL).tex

$(MANUAL).bbl: $(MANUAL).bib
	$(LATEX) $(MANUAL).tex
	$(BIBTEX) $(MANUAL)

$(MANUAL).bib: $(MANUAL).bst $(ABBREV)


clean:
	$(RM) *.aux *.log *.rel *.bbl *.blg;\
	$(RM) -r .auto/

distclean: clean
	$(RM) *.dvi *.ps *.pdf
