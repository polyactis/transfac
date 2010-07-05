#!/usr/bin/env python
"""
Usage: ConstructTransfacProf.py [OPTIONS] -p embl_profile -m matrix_file -o output_file

Option:
	-p ...	EMBL-format profile
	-m ...	EMBL-format matrix data
	-o ...	output file
	--mc=...	matrix similarity score cutoff
	--cc=...	core similarity score cutoff
	-b, --debug	DEBUG
	-r, --report	report the progress (time before each query)
	-h, --help	show this help

Examples:

Description:
	This program transforms the embl_profile into the match-format profile.
	matrix_file is needed to find the matrix accession.
	If mc and cc are not given, just use the corresponding part in embl_profile.
"""

import sys, getopt, os, cStringIO
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
sys.path += [os.path.join(os.path.expanduser('~/script/microarray/bin'))]
from MdbId2UnigeneId import unigene_data_block_iterator

class ConstructTransfacProf:
	"""
	09-10-05
	"""
	def __init__(self, embl_profile, matrix_file, output_file, mc, cc, debug=0, report=0):
		"""
		09-10-05
		"""
		self.embl_profile = embl_profile
		self.matrix_file = matrix_file
		self.output_file = output_file
		self.mc = mc
		self.cc = cc
		self.debug = int(debug)
		self.report = int(report)
	
	def get_matrix_id2acc(self, matrix_file):
		"""
		09-10-05
		
		"""
		sys.stderr.write("Setting up matrix_id2acc from %s..."%matrix_file)
		inf = open(matrix_file,'r')
		iter = unigene_data_block_iterator(inf)
		matrix_id2acc = {}
		for block in iter:
			if block=='':	#the last nothing block
				break
			block = cStringIO.StringIO(block)
			acc = None
			id = None
			for line in block:
				if line.find('AC'+' '*2)==0:
					acc = line[4:-1]
				if line.find('ID'+' '*2)==0:
					id = line[4:-1]
			if acc and id:	#the first block of the matrix_file has no matrix
				matrix_id2acc[id] = acc
		
		del iter, inf
		sys.stderr.write("Done\n")
		return matrix_id2acc
	
	def parse_embl_profile(self, embl_profile, output_file, matrix_id2acc, mc_given, cc_given):
		"""
		09-10-05
		09-11-05
			the firstline of output_file changed to the basename of the output_file, not the full path.
			The format is based on internal profile(minSUM92.prf), which is a little bit different from
				what is said in the documentation file from Kangyu(No first blank).
		"""
		sys.stderr.write("Parsing from %s to %s..."%(embl_profile, output_file))
		inf = open(embl_profile,'r')
		iter = unigene_data_block_iterator(inf)
		of = open(output_file, 'w')
		#write some header information
		of.write("%s\n"%os.path.basename(output_file))
		of.write("From %s.\n"%(os.path.basename(embl_profile)))
		of.write(" MIN_LENGTH 300\n")
		of.write("0.0\n")
		for block in iter:
			if block=='':	#the last nothing block
				break
			block = cStringIO.StringIO(block)
			id = None
			mc = None
			cc = None
			for line in block:
				if line.find('ID'+' '*1)==0:
					id = line[3:-1]
				if line.find('MC'+' '*1)==0:
					mc = line[3:-1]
				if line.find('CC'+' '*1)==0:
					cc = line[3:-1]
			if mc_given:
				mc = mc_given
			if cc_given:
				cc = cc_given
			if id and mc and cc:
				acc = matrix_id2acc[id]
				of.write(' 1.000000 %s %s %s %s\n'%(cc,mc, acc,id))	#first character is blank
		of.write('//\n')
		del inf
		sys.stderr.write("Done\n")
	
	def run(self):
		"""
		09-10-05
		"""
		matrix_id2acc = self.get_matrix_id2acc(self.matrix_file)
		self.parse_embl_profile(self.embl_profile, self.output_file, matrix_id2acc, self.mc, self.cc)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["mc=", "cc=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "p:m:o:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	embl_profile = None
	matrix_file = None
	output_file = None
	mc = None
	cc = None
	debug = 0
	report = 0
	commit = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-p",):
			embl_profile = arg
		elif opt in ("-m",):
			matrix_file = arg
		elif opt in ("-o",):
			output_file = arg
		elif opt in ("--mc",):
			mc = float(arg)
		elif opt in ("--cc",):
			cc = float(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if embl_profile and matrix_file and output_file:
		instance = ConstructTransfacProf(embl_profile, matrix_file, output_file, mc, cc, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
