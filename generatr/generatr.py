##
## Generic imports
import os
import sys
import argparse
import logging as log
import pkg_resources
from itertools import product

##
## Subpackage/s??
from .dtdvalidate.validation import Colour as clr
from .dtdvalidate.validation import ConfigReader

class Generatr:
	def __init__(self):

		"""
		Simple data scraping script which takes a sam file as input
		and collects the top 3 references in regards to aligned read count.
		"""
		##
		## Package data
		self.package_exampleXML = pkg_resources.resource_filename(__name__, 'dtdvalidate/example_input.xml')
		self.package_configDTD = pkg_resources.resource_filename(__name__, 'dtdvalidate/xml_rules.dtd')


		##
		## Argument parser from CLI
		self.parser = argparse.ArgumentParser(prog='generatr',description='RefGeneratr: Dynamic multi-loci/mutli-repeat tract microsatellite sequence generator.')
		self.parser.add_argument('-i','--input',help='Input data. Path to input XML document with desired sequence information.',nargs=1,required=True)
		self.parser.add_argument('-o','--output',help='Output path. Specify a directory wherein your *.fa reference will be saved.',nargs=1,required=True)
		self.parser.add_argument('-v','--verbose',help='Verbose mode. Enables terminal user feedback.',action='store_true')
		self.args = self.parser.parse_args()

		##
		## Sets up input directories and verbose mode if requested
		self.input_directory = self.args.input[0]
		self.output_directory = self.args.output[0]
		if self.args.verbose:
			log.basicConfig(level=log.DEBUG, format="%(message)s")

		##
		## User feedback and run application
		## Checks for appropriate input/output -- if self.iocheck() returns false, quit
		log.info('{}{}{}{}'.format(clr.bold,'gtr__ ',clr.end,'RefGeneratr: microsatellite reference sequence generator.'))
		log.info('{}{}{}{}'.format(clr.bold,'gtr__ ',clr.end,'alastair.maxwell@glasgow.ac.uk'))
		if not self.iocheck():
			log.error('{}{}{}{}'.format(clr.red,'gtr__ ',clr.end,'Exiting..'))
			sys.exit(2)

		##
		## Each input sam will have an entry in a results dictionary which is then output to
		## the master output folder in a single CSV file for easy interpretation
		self.greeting = 'hi'

	def iocheck(self):
		"""
		Checks input isfile/isXML;;	Checks exampleXML/DTD validates
		Checks inputXML validates;; Checks output file
		Return True if all OK, else False
		"""

		##
		## Check that specified input is an XML document
		if not (os.path.isfile(self.input_directory)):
			log.error('{}{}{}{}'.format(clr.red,'gtr__ ',clr.end,'I/O: Specified input path is not a file!'))
			return False
		if not (self.input_directory.endswith('.xml') or self.input_directory.endswith('.XML')):
			log.error('{}{}{}{}'.format(clr.red,'gtr__ ',clr.end,'I/O: Specified input file is not XML!'))
			return False

		##
		## Validate exampleXML against package DTD
		ConfigReader(self.package_configDTD, self.package_exampleXML)

		##
		## Validate input XML against package DTD

		##
		## Specified output check
		if not (self.output_directory.endswith('.fasta') or self.output_directory.endswith('.fa') or self.output_directory.endswith('.fas')):
			log.error('{}{}{}{}'.format(clr.red,'gtr__ ',clr.end,'I/O: Specified output path is not targetting a FASTA file!'))
			return False

		return True

	def workflow(self):
		print self.greeting


		default_locus = '>Test_locus_'
		fiveprime_flank = 'fiveprimeflank'
		threeprime_flank = 'threeprimeflank'

		input_repeat_regions = [
			{'label': 'CTG',
			 'start': 1,
			 'end': 20},

			{'label': 'CTG',
			 'start': 1,
			 'end': 30},

			{'label': 'CTG',
			 'start': 1,
			 'end': 20}
		]

		intervening_regions = [
			{'first': 'CTG',
			 'intervening': 'CCGCCG'},

			{'first': 'CTG',
			 'intervening': 'CCG'}

		]

		range_lists = []

		# Process the range for each of these repeat regions
		for region in input_repeat_regions:
			region['range'] = range(region['start'], region['end']+1)
			range_lists.append(region['range'])

		# For each list of label ranges, generate the cartesian product.
		generator = product(*range_lists)

		# This is going to be a sorted list, we start sorting from the first element of each tuple.
		s = list(product(*range_lists))

		# Do the sorting at each position. This is probably inefficient, and there's likely a better way of doing this...
		for i in range(1, len(input_repeat_regions)):
			s = sorted(s, key=lambda x: x[i])

		# Make the strings.

		def get_sandwich(first_region):
			"""
			If a intervening sequence exists just after the first_region, returns that intervening sequence as a string.
			If no intervening sequence exists, an emptystring is returned.
			"""
			for i in intervening_regions:
				if i['first'] == first_region:
					return i['intervening']

			return ''

		# Do the logic

		for range_tuple in s:  # range_tuple is the number of times each repeat regions must be dumped.
			i = 0
			big_string = ''

			locus_label = ''
			for num_times in range_tuple:
				label = input_repeat_regions[i]['label']
				sandwich = get_sandwich(label)

				repeat_region = label * num_times
				complete_region = '{repeat_region}{sandwich}'.format(repeat_region=repeat_region,sandwich=sandwich)
				big_string = '{big_string}{complete_region}'.format(big_string=big_string,complete_region=complete_region)

				#print label, num_times,
				locus_label += label + '_' + str(num_times) + '_'



				i += 1

			print

			entire_label = default_locus + locus_label
			entire_label = entire_label[:-1]
			print entire_label
			print '{fiveprime_flank}{big_string}{threeprime_flank}'.format(fiveprime_flank=fiveprime_flank,
																		   big_string=big_string,
																		   threeprime_flank=threeprime_flank)

			locus_label = ''



def main():
	Generatr()