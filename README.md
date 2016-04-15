RefGeneratr: Dynamic multi-loci/mutli-repeat tract microsatellite sequence generator
=========================================================
SBMLCartographer is a python package which utilises SBML API's to track the transition
of Carbon-13 atoms between biomass reactions in order to analysis metabolic flux (typically in Yeast).
This was started as a MSc level thesis project, but is being re-written and distributed as a package.
It is in very early stages and has no functionality as of present. Constant work in progress.

RefGeneratr (generatr) is a python script/package which generates a reference genetic sequence (*.fasta) for use in sequence alignment.
Microsatellite repeat regions can vary in scope and loci count, so this software has the ability to dynamically handle an undetermined
amount of repeat regions within each loci, with intervening sequences if desired. Endusers can specify as many regions/loci as desired, through
a simple XML document. This is parsed, and output in the standard *.fasta format is provided.

There should be no non-standard requirements, if your operating system has Python 2.7 installed. Setuptools/PIP should handle any missing python
packages, however.

What's New
==========
Pre-release status.
Tidying up algorithm into distributable.


Installation Prerequisites
==========================
Work in progress, so don't count on this being filled out any time soon.
Should just be..
    $ python setup.py install
or
    $ pip install RefGeneratr
..when development is complete.


Usage
=====

Here's how generatr will be controlled when development is complete:

    $ generatr [-v/--verbose] [-i/--input <Path to input.xml>] [-o/--output <Desired *.fasta file output>]

Should be fairly self-explanatory.











