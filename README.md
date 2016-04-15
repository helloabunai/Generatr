RefGeneratr: Dynamic multi-loci/multi-repeat tract microsatellite reference sequence generator
=========================================================
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











