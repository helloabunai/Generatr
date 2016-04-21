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
		log.info('{}{}{}{}{}'.format('\n', clr.bold,'gtr__ ',clr.end,'RefGeneratr: microsatellite reference sequence generator.'))
		log.info('{}{}{}{}'.format(clr.bold,'gtr__ ',clr.end,'alastair.maxwell@glasgow.ac.uk'))
		self.input_dictionary = {}
		if not self.iocheck():
			log.error('{}{}{}{}'.format(clr.red,'gtr__ ',clr.end,'Exiting..'))
			sys.exit(2)

		##
		## Loop over every loci that we scraped from XML
		## Extract info and format tailor to generator methods
		## Save reference string to structure, when all complete, pass to file writer
		self.loci_strings = []
		for k,v in self.input_dictionary.iteritems():

			##
			## Single Loci mode, structure returned an individual dictionary
			if type(v) == dict:
				log.info('{}{}{}{}'.format(clr.bold,'gtr__ ',clr.end,'Detected a single locus. Processing..'))
				locus_dictionary = self.loci_collector(v)
				locus_string = self.generate_loci_reference(locus_dictionary)
				self.loci_strings.append(locus_string)
				self.write_output()
			##
			## Multi loci mode, structure returned a list of dictionaries
			if type(v) == list:
				log.info('{}{}{}{}'.format(clr.bold,'gtr__ ',clr.end,'Detected multiple loci. Processing..'))
				for raw_locus in v:
					locus_dictionary = self.loci_collector(raw_locus)
					locus_string = self.generate_loci_reference(locus_dictionary)
					self.loci_strings.append(locus_string)
				self.write_output()

		log.info('bye')

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
		## Validate exampleXML, then inputXML against package DTD
		ConfigReader(self.package_configDTD, self.package_exampleXML)
		self.input_dictionary = ConfigReader(self.package_configDTD, self.input_directory).return_dict()

		##
		## Specified output check
		if not (self.output_directory.endswith('.fasta') or self.output_directory.endswith('.fa') or self.output_directory.endswith('.fas')):
			log.error('{}{}{}{}'.format(clr.red,'gtr__ ',clr.end,'I/O: Specified output path is not targetting a FASTA file!'))
			return False

		return True

	@staticmethod
	def loci_collector(raw_locus):

		"""
		Function to extract (now verified) information to pass to cartesian string generator
		Assuming all information is now valid so no more checks are carried out
		Not the most efficient way of doing it.. but hey
		"""

		##
		## FastA dictates each reference should begin with '>'
		## So check and adhere if required
		loci_label = raw_locus['@label']
		if not loci_label.startswith('>'):
			loci_label = '>' + loci_label

		##
		## Flanks to be filled in
		fiveprime_flank = ''
		threeprime_flank = ''

		##
		## Repeat region(s) to be filled in
		input_repeat_regions = []

		##
		## Intervening sequence(s) to be filled in
		intervening_regions = []

		##
		## Dictionaries for each reference region
		loci_data = raw_locus['input']
		for sequence_parameters in loci_data:

			##
			## Scraping 5' flank
			if sequence_parameters['@type'] == 'fiveprime':
				fiveprime_flank = sequence_parameters['@flank']

			##
			## Scraping 3' flank
			if sequence_parameters['@type'] == 'threeprime':
				threeprime_flank = sequence_parameters['@flank']

			##
			## Scraping repeat regions
			if sequence_parameters['@type'] == 'repeat_region':
				current_region = {'unit':sequence_parameters['@unit'],
								  'start':sequence_parameters['@start'],
								  'end':sequence_parameters['@end'],
								  'order':sequence_parameters['@order']}
				input_repeat_regions.append(current_region)

			##
			## Scraping intervening regions
			if sequence_parameters['@type'] == 'intervening':
				current_interv = {'prior':sequence_parameters['@prior'],
								  'sequence':sequence_parameters['@sequence']}
				intervening_regions.append(current_interv)

		locus_dictionary = {'label': loci_label,
						   '5P_flank': fiveprime_flank,
						   'repeat_regions': input_repeat_regions,
						   'intervening_regions': intervening_regions,
						   '3P_flank': threeprime_flank}

		return locus_dictionary

	@staticmethod
	def generate_loci_reference(locus_dictionary):

		##
		## Parameters for the current locus
		## that we got passed in locus_dictionary
		loci_label = locus_dictionary['label']
		five_prime = locus_dictionary['5P_flank']
		three_prime = locus_dictionary['3P_flank']
		repeat_regions = locus_dictionary['repeat_regions']
		intervening_regions = locus_dictionary['intervening_regions']

		##
		## Process the explicit ranges for each repeat region
		possible_ranges = []
		for region in repeat_regions:
			region['range'] = range(int(region['start']), int(region['end'])+1)
			possible_ranges.append(region['range'])

		##
		## For each list of label ranges, generate cartesian product
		## Sort the list into our desired order..
		## start based on first tuple element
		## Then, do sorting at each position (inefficient hhhhhehe)
		s = list(product(*possible_ranges))
		for i in range(1, len(repeat_regions)):
			s = sorted(s, key=lambda x:x[i])

		## Now we make the strings
		## i.e. joining intervening sequence to their preceding repeat region
		## assuming that intervening's "prior" flag == specified order flag
		def wedge_intervening(region_order):
			for section in intervening_regions:
				if section['prior'] == region_order:
					return section['sequence']
			return ''

		##
		## Generation of entire reference string
		complete_reference_string = ''
		for range_tuple in s:
			i = 0
			nonflank_sequences = ''
			haplotype_label = ''

			##
			## For this specific range tuple, how many times do we want the repeat unit to occur?
			## organise this string with flanks, complete string for this individual 'reference point'
			for num_times in range_tuple:

				##
				## Information on reference characteristics
				repeat_unit = repeat_regions[i]['unit']
				repeat_order = repeat_regions[i]['order']
				post_repeatregion = wedge_intervening(repeat_order)

				##
				## Generating the string as directed
				repeat_region = repeat_unit * num_times
				complete_region = '{}{}'.format(repeat_region,post_repeatregion)
				nonflank_sequences = '{}{}'.format(nonflank_sequences,complete_region)

				##
				## A label anchor appended to the given label
				## conveying the specified information for the current reference
				haplotype_label += '_' + repeat_unit + str(num_times)
				i += 1

			##
			## Combine data, append to entire locus' string
			reference_label = loci_label + haplotype_label
			reference_sequence = '{}{}{}'.format(five_prime, nonflank_sequences, three_prime)
			reference_string = '{}\n{}'.format(reference_label, reference_sequence)
			complete_reference_string += reference_string

		return complete_reference_string

	def write_output(self):

		print len(self.loci_strings)

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