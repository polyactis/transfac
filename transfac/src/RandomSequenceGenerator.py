#!/usr/bin/env python
"""
Usage: RandomSequenceGenerator.py -o output_dir [OPTIONS]

Option:
	-o ...,	output_dir
	-a ...,	A percentage, 0.25(default)
	-t ...,	T percentage, 0.25(default)
	-g ...,	G percentage, 0.25(default)
	-n ...,	no_of_files, 10(default)
	-p ...,	output_prefix, 'random_sequence'(default)
	-s ...,	no_of_sequences per file, 1000(default)
		TRANSFAC's match has upper limit 2G output file size
	-l ...,	length of each sequence, 10000(default)
	-y ...,	(IGNORE)
	-b	enable debugging, no debug by default
	-r	report the progress(a number)
	-h, --help              show this help

Examples:
	RandomSequenceGenerator.py -o /tmp/test -a 0.3 -t 0.1 -g 0.2

Description:
	Generate random sequences in fasta format for TRANSFAC.
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
import random, getopt

class RandomSequenceGenerator:
	def __init__(self, output_dir=None, A_percentage=0.25, T_percentage=0.25, G_percentage=0.25,\
		no_of_files=10000, output_prefix='random_sequence', no_of_sequences_per_file=1000, \
		length_of_sequence=10000, type=1, debug=0, report=0):
		self.output_dir = output_dir
		self.A_percentage = float(A_percentage)
		self.T_percentage = float(T_percentage)
		self.G_percentage = float(G_percentage)
		self.no_of_files = int(no_of_files)
		self.output_prefix = output_prefix
		self.no_of_sequences_per_file = int(no_of_sequences_per_file)
		self.length_of_sequence = int(length_of_sequence)
		self.type = int(type)
		self.debug = int(debug)
		self.report = int(report)
	
	def generate_random_segment(self, base_choice, length):
		"""
		01-16-06
			based on base_choice(either AT or GC), generate a segment
			50% possibility to choose from base_choice
		01-19-06
			defunct
		"""
		random_segment = []
		for i in range(length):
			random_segment.append(base_choice[random.randint(0,1)])
		return random_segment
	
	def generate_one_file(self, output_fname, no_of_sequences_per_file, A_percentage, T_percentage, G_percentage, length_of_sequence):
		"""
		01-16-06
			fix a bug in output format
		01-19-06
			add A_percentage, T_percentage, G_percentage
		"""
		sys.stderr.write("Generating for %s..."%os.path.basename(output_fname))
		outf = open(output_fname, 'w')
		A_segment_length = int(A_percentage*length_of_sequence)
		T_segment_length = int(T_percentage*length_of_sequence)
		G_segment_length = int(G_percentage*length_of_sequence)
		C_segment_length = length_of_sequence - A_segment_length - T_segment_length - G_segment_length
		total_segment = ['A']*A_segment_length
		total_segment += ['T']*T_segment_length
		total_segment += ['G']*G_segment_length
		total_segment += ['C']*C_segment_length
		
		for i in range(no_of_sequences_per_file):
			self.counter += 1
			"""
			GC_segment = self.generate_random_segment(['G', 'C'], GC_segment_length)
			AT_segment = self.generate_random_segment(['A', 'T'], length_of_sequence-GC_segment_length)
			total_segment = GC_segment + AT_segment
			"""
			random.shuffle(total_segment)	#random shuffle the AT and GC content
			outf.write('>%s\n'%self.counter)
			outf.write("%s\n"%''.join(total_segment))	#01-16-06
		del outf
		sys.stderr.write("Done.\n")
	
	def run(self):
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		self.counter = 0	#as sequence id
		for i in range(self.no_of_files):
			output_fname = os.path.join(self.output_dir, '%s%s'%(self.output_prefix, i))
			self.generate_one_file(output_fname, self.no_of_sequences_per_file, self.A_percentage, self.T_percentage, \
				self.G_percentage, self.length_of_sequence)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "debug"]
	opts, args = getopt.getopt(sys.argv[1:], "ho:a:t:g:n:p:s:l:y:br", long_options_list)
	
	output_dir = None
	A_percentage = 0.25
	T_percentage = 0.25
	G_percentage = 0.25
	no_of_files = 10
	output_prefix = 'random_sequence'
	no_of_sequences_per_file = 1000
	length_of_sequence = 10000
	type = 1
	debug = 0
	report = 0

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-o",):
			output_dir = arg
		elif opt in ("-a",):
			A_percentage = float(arg)
		elif opt in ("-t",):
			T_percentage = float(arg)
		elif opt in ("-g",):
			G_percentage = float(arg)
		elif opt in ("-n",):
			no_of_files = int(arg)
		elif opt in ("-p",):
			output_prefix = arg
		elif opt in ("-s",):
			no_of_sequences_per_file = int(arg)
		elif opt in ("-l",):
			length_of_sequence = int(arg)
		elif opt in ("-y",):
			type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r",):
			report = 1
		
	if output_dir and G_percentage:
		instance = RandomSequenceGenerator(output_dir, A_percentage, T_percentage, G_percentage, no_of_files, output_prefix, \
			no_of_sequences_per_file, length_of_sequence, type, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
