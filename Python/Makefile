
.PHONY: help all build install doc doc-html doc-pdf
help:
	@echo ""
	@echo "Please use 'make <target>' where <target> is one of"
	@echo ""
	@echo "   all          to build this package and run the testsuite"
	@echo "   build        to build this package"
	@echo "   install      to install this package"
	@echo "   test         to run the unit-test suite"
	@echo "   doc"
	@echo "   doc-html     to build HTML documentation"
	@echo "   doc-pdf      to build the PDF codumentation"
	@echo ""

all: build test

build:
	python setup.py build
	python setup.py build_ext --inplace

install: build
	python setup.py install

clean:
	python setup.py clean
	rm -Rf build
	rm -Rf doc/.build

test: build
	nosetests

doc: doc-html

doc-html:
	cd doc; make html

doc-pdf:
	cd doc; make latex
	cd doc/.build/latex; make all-pdf

