import sys
import logging as log
from collections import defaultdict
from xml.etree import cElementTree
from lxml import etree

class Colour:

	def __init__(self):
		pass

	purple = '\033[95m'
	cyan = '\033[96m'
	darkcyan = '\033[36m'
	blue = '\033[94m'
	green = '\033[92m'
	yellow = '\033[93m'
	red = '\033[91m'
	bold = '\033[1m'
	underline = '\033[4m'
	end = '\033[0m'

class ConfigReader(object):

	"""
	The configuration file reader.
	Opens a configuration file, and if valid, converts the parameters within the file to a dictionary object,
	reader to be viewed through accessing the config_dict variable.
	"""

	def __init__(self, dtdfile, config_filename=None):

		##
		## Instance variables
		self.config_filename = config_filename
		self.dtd_filename = dtdfile

		##
		## Check for configuration file (just incase)
		if self.config_filename is None:
			log.error("No configuration file specified!")
		else:
			self.config_file = etree.parse(self.config_filename)

		##
		## Check config vs dtd, parse info to dictionary, validate vs ruleset
		self.validate_against_dtd()
		self.set_dictionary()
		self.trigger = self.validate_config()

	def validate_against_dtd(self):

		"""
		Validate input config against DTD ruleset
		i.e. confirms conformation of XML structure
		"""

		##
		## Open > etree.DTD object
		dtd_file = open(self.dtd_filename, 'r')
		dtd_object = etree.DTD(dtd_file)

		##
		## If validation fails, close the object (memory) and raise an error
		if not dtd_object.validate(self.config_file):
			dtd_file.close()
			log.error("DTD validation failure {0}: {1}".format(self.config_filename, dtd_object.error_log.filter_from_errors()[0]))
			sys.exit(2)
		dtd_file.close()

	def set_dictionary(self):

		"""
		Takes the now validated XML and extracts information from the tree into
		a python dictionary {key: value}. This dictionary will be used for variables
		within the pipeline. Recursion adapted from http://stackoverflow.com/a/9286702
		"""

		def recursive_generation(t):

			d = {t.tag: {} if t.attrib else None}
			children = list(t)

			##
			## If list was populated, create dictionary, Append keys
			if children:
				dd = defaultdict(list)

				for dc in map(recursive_generation, children):
					for k, v in dc.iteritems():
						dd[k].append(v)
				d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.iteritems()}}

			##
			## Values for key
			if t.attrib:
				d[t.tag].update(('@' + k, v) for k, v in t.attrib.iteritems())

			if t.text:
				text = t.text.strip()
				if children or t.attrib:
					if text:
						d[t.tag]['#text'] = text
				else:
					d[t.tag] = text
			return d

		##
		## Takes the formatted xml doc, puts through generator, returns dictionary
		string_repr = etree.tostring(self.config_file, pretty_print=True)
		element_tree = cElementTree.XML(string_repr)

		self.config_dict = recursive_generation(element_tree)
		self.config_dict = self.config_dict[self.config_dict.keys()[0]]

	def validate_config(self):

		"""
		Method which validates the configuration file's contents.
		If all pass, guarantees that the settings dictionary is full of valid settings!
		"""
		valid_types = ['fiveprime','repeat_region','intervening','threeprime']
		valid_bases = ['A','T','G','C','U','N']

		def parameter_checker(input_data):

			subtrigger = False

			for sequence_parameters in input_data:

				## Check type integrity
				param_type = sequence_parameters['@type']
				if not param_type in valid_types:
					log.error('{}{}{}{}'.format(Colour.red,'gtr__ ',Colour.end,'CFG: Invalid parameter type detected in XML config.'))

				##
				## Test five prime flank integrity
				if param_type == 'fiveprime':
					fiveprime_flank = sequence_parameters['@flank']
					for character in fiveprime_flank:
						if not character in valid_bases:
							log.error('{}{}{}{}'.format(Colour.red,'gtr__ ',Colour.end,'CFG: Invalid base character detected in five prime flank.'))
							subtrigger=True

				##
				## Test repeat region(s) integrity
				if param_type == 'repeat_region':
					repeat_region = sequence_parameters['@unit']
					for character in repeat_region:
						if not character in valid_bases:
							log.error('{}{}{}{}'.format(Colour.red,'gtr__ ',Colour.end,'CFG: Invalid base character detected in repeat unit.'))
							subtrigger=True

					region_start = sequence_parameters['@start']
					region_end = sequence_parameters['@end']
					region_order = sequence_parameters['@order']
					for region in [region_start, region_end, region_order]:
						if not region.isdigit():
							log.error('{}{}{}{}'.format(Colour.red,'gtr__ ',Colour.end,'CFG: Non-integer character detected in repeat region start/end/order.'))
							subtrigger=True

				##
				## Test intervening sequence(s) integrity
				if param_type == 'intervening':
					intevening_seq = sequence_parameters['@sequence']
					for character in intevening_seq:
						if not character in valid_bases:
							log.error('{}{}{}{}'.format(Colour.red,'gtr__ ',Colour.end,'CFG: Invalid base character found in intervening sequence.'))
							subtrigger=True

				##
				## Test three prime flank integrity
				if param_type == 'threeprime':
					threeprime_flank = sequence_parameters['@flank']
					for character in threeprime_flank:
						if not character in valid_bases:
							log.error('{}{}{}{}'.format(Colour.red,'gtr__ ',Colour.end,'CFG: Invalid base character detected in three prime flank.'))
							subtrigger=True

			return subtrigger

		trigger = False
		for k,v in self.config_dict.iteritems():

			if type(v) == dict:
				loci_data = v['input']
				trigger = parameter_checker(loci_data)

			if type(v) == list:
				for loci_entry in v:
					loci_data = loci_entry['input']
					trigger = parameter_checker(loci_data)

		if trigger:
			log.error('{}{}{}{}'.format(Colour.red,'gtr__ ',Colour.end,'CFG: XML parameter validation failure. Exiting.'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(Colour.green,'gtr__ ',Colour.end,'CFG: XML Validation success!'))

	def return_dict(self):
		return self.config_dict