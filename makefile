SHELL := /bin/bash

.PHONY: all packages main

all: packages main compile

packages:
	@julia code/packages.jl

main:
	@julia code/main.jl --verbose

compile:
	TEXFILE := figures_and_tables
	cd ./draft; \
	pdflatex $(TEXFILE).tex; \
	pdflatex $(TEXFILE).tex; \
	rm -f $(TEXFILE).aux $(TEXFILE).bbl $(TEXFILE).blg $(TEXFILE).log $(TEXFILE).out $(TEXFILE).toc

