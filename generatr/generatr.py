from itertools import product

class Generatr:
	def __init__(self):
		self.greeting = 'hi'

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