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
		self.validate_config()

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

		##TODO generate proper loci-based dictionary (one key atm???)

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
		trigger = False

		for k, v in self.config_dict.iteritems():
			print 'KEY: ', k
			for thing in v:
				print 'LOCI LABEL: ', thing['@label']
				print 'LOCI DATA: ', thing['input']
				print '\n'

		##TODO validate input


		###
		### Main configuration instance settings
		#
		#data_directory = self.config_dict['@data_dir']
		#if not os.path.exists(data_directory):
		#	log.error('{}{}{}{}'.format(Colour.red, 'shd__', Colour.end, 'XML Config: Specified data directory could not be found.'))
		#	trigger = True
		#for fqfile in glob.glob(os.path.join(data_directory, '*')):
		#	if not (fqfile.endswith('.fq') or fqfile.endswith('.fastq') or fqfile.endswith('.fq.gz') or fqfile.endswith('.fastq.gz')):
		#		log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Non FastQ/GZ data detected in specified input directory.'))
		#		trigger = True
		#reference_directory = self.config_dict['@reference_file']
		#if not os.path.isfile(reference_directory):
		#	log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified reference file could not be found.'))
		#	trigger = True
		#if not (reference_directory.endswith('.fa') or reference_directory.endswith('.fas')):
		#	log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified reference file is not an fa/fas file.'))
		#	trigger = True
		#
		###
		### Instance flag settings
		#
		#sequence_qc_flag = self.config_dict['instance_flags']['@quality_control']
		#if not (sequence_qc_flag == 'True' or sequence_qc_flag == 'False'):
		#	log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Sequence Quality control flag is not set to True/False.'))
		#	trigger = True
		#alignment_flag = self.config_dict['instance_flags']['@sequence_alignment']
		#if not (alignment_flag == 'True' or alignment_flag == 'False'):
		#	log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Sequence Alignment flag is not set to True/False.'))
		#	trigger = True
		#genotype_flag = self.config_dict['instance_flags']['@genotype_prediction']
		#if not (genotype_flag == 'True' or genotype_flag == 'False'):
		#	log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Genotype Prediction control flag is not True/False.'))
		#	trigger = True
		#
		#if sequence_qc_flag == 'True' and alignment_flag == 'False':
		#	log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Quality control selected, but not alignment. Invalid selection.'))
		#	trigger = True
		#
		#if trigger:
		#	log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Failure, exiting.'))
		#	sys.exit(2)
		#else:
		#	log.info('{}{}{}{}'.format(Colour.green, 'shd__ ', Colour.end, 'XML Config: Parsing parameters successful!'))
