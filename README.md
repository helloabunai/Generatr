RefGeneratr: Dynamic multi-loci/multi-repeat tract microsatellite reference sequence generator
==============================================================================================
RefGeneratr (generatr) is a python script/package which generates a reference genetic sequence (*.fasta) for use in sequence alignment.
Microsatellite repeat regions can vary in scope and loci count, so this software has the ability to dynamically handle an undetermined
amount of repeat regions within each loci, with intervening sequences if desired. Endusers can specify as many regions/loci as desired, through
a simple XML document. This is parsed, and output in the standard *.fasta format is provided.

Generatr requires lxml, which setuptools should install for you during setup.

What's New
==========
Everything


Installation Prerequisites
==========================

Assuming that lxml is installed, or you wish setuptools to handle installation for you, the following should suffice. For now, download the source and run:

    $ python setup.py install

You may or may not required sudo, it depends on your system. This will install the package for you, so it can be launched with 'generatr' from the command line.
Eventually, the package will be uploaded onto PIP so that you can install directly from a terminal.

Hardware Requirements
=====================

Nothing spectacular, any computer should run it fine. However, if you desire to generate a reference with a large amount of repeat units and/or loci, available
system memory will be a bottleneck.

Usage
=====

Here's how to use generatr:

    $ generatr [-v/--verbose] [-i/--input <Path to input.xml>] [-o/--output <Desired *.fasta file output>]

-v enables terminal user feedback.

-i is a path to an XML file containing your desired information, which adheres to the requirements outlined below.

-o is a path to your desired output *.fasta/*.fa/*.fas file.

XML Requirements
=====

An example XML file is as follows:

    <?xml version="1.0"?>
    <data>
        <loci label="example_loci_one">
            <input type="fiveprime" flank="GCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC"/>
            <input type="repeat_region" order="1" unit="CAG" start="1" end="100"/>
            <input type="intervening" sequence="CAACAGCCGCCA" prior="1"/>
            <input type="repeat_region" order="2" unit="CCG" start="1" end="20"/>
            <input type="threeprime" flank="CCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCT"/>
        </loci>
    </data>

The input regions have been made as straight forward as possible. If you desire multiple loci within one reference file,
additional <loci> tags should be presented, with the respective sequence parameters nested within. There is technically no limitation
on how many loci you can specify, although testing has not gone beyond any reasonable figures.

The possible sequence parameters are as follows:

    <input type="fiveprime" flank="<string>"/>

This is the input for a five prime flank sequence. The 'type' must be 'fiveprime', and any valid sequence can be present within
the 'flank' variable. Valid sequence is a string that consists of A,G,C,T,U,N. No other characters are considered valid.

    <input type="repeat_region" order="<integer>" unit="<string>" start="<integer>" end="<integer>"/>

This is the input for a repeat region. The order flag indicates where in the 'sequence' it resides. Unit equates to the repeated unit
of sequence, and start/end are integers for the range you wish this repeat unit to repeat over. Generatr is useful as it can handle an unspecified
number of repeat regions for each loci.

    <input type="intervening" sequence="<string>" prior="<integer>"

The intervening flag is for interrupted repeat sequences. Your intervening sequence is specified under 'sequence', and the repeat_region
which this intervening sequence follows, is indicated in 'prior'. E.G. if an intervening sequence follows a repeat_region that was order="1",
the intervening prior value would also be "1". Generatr can handle zero, one or multiple intervening sequences; the only stipulation for the sequence
to appear correctly is for the user to accurately input the preceeding repeat_region's 'order' value under the respective intervening region's 'prior' value.

    <input type="threeprime" flank="<string>"/>

The input for a three prime flank follows the same logic as described for five prime.

Thanks for reading. If you have any questions or trouble with installation, please feel free to e-mail me at alastair[dot]maxwell[at]glasgow[dot]ac[dot]uk.












