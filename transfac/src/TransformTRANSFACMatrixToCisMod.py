#!/usr/bin/env python
"""
Usage: TransformTRANSFACMatrixToCisMod.py -p PROFILE_FNAME -m MATRIX_FNAME -o OUTPUT_FNAME [OPTIONS]

Option:
	-p ...,	profile(from TRANSFAC) filename
	-m ...,	matrix.dat filename from TRANSFAC
	-o ...,	output filename
	-b,	enable debug
	-r,	enable report
	-h, --help	show this help

Examples:
	TransformTRANSFACMatrixToCisMod.py -p data/prfs/vertebrates_minSUM_highQual.prf -m data/matrix.dat -o /tmp/vertebrates_minSUM_highQual.pwm
	
Description:
	Program to convert TRANSFAC's matrix.dat to CisModScan's pwm format.
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.expanduser('~/script64/microarray/bin'))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.expanduser('~/script/microarray/bin'))
import getopt, re, cStringIO
from sets import Set
from MdbId2UnigeneId import unigene_data_block_iterator

class TransformTRANSFACMatrixToCisMod:
	def __init__(self, profile_fname=None, matrix_fname=None, output_fname=None, debug=0, report=0):
		"""
		2006-08-11
		"""
		self.profile_fname = profile_fname
		self.matrix_fname = matrix_fname
		self.output_fname = output_fname
		self.debug = int(debug)
		self.report = int(report)
	
	def getProfileIdSet(self, profile_fname):
		sys.stderr.write("Getting a profile set...")
		profile_id_set = Set()
		inf = open(profile_fname)
		for line in inf:
			if line[:2]=='ID':
				profile_id_set.add(line[3:-1])
		del inf
		sys.stderr.write("Done.\n")
		return profile_id_set
	
	def transformMatrix(self, matrix_fname, output_fname, profile_id_set):
		"""
		2006-08-14
			similar to transfacdb.py's transfac_embl_parse
		"""
		sys.stderr.write("Transforming matrix...")
		inf = open(matrix_fname,'r')
		iter = unigene_data_block_iterator(inf)
		outf = open(output_fname, 'w')
		profile_id_outf = open('%s.id_mapping'%output_fname, 'w')
		block_no = 0
		pwm_line_pattern = re.compile(r'\d\d  ')	#used to identify pwm lines, first two characters are numbers
		profile_id_list = []
		divid_f= lambda x: x/sum(pwm_row)
		for block in iter:
			block_no += 1
			if block_no == 1:	#skip the first block
				continue
			pwm_list = []
			block = cStringIO.StringIO(block)
			for line in block:
				try:
					if pwm_line_pattern.match(line):
						ls = line.split()
						pwm_row = map(float, ls[1:-1])
						pwm_row = map(divid_f, pwm_row)
						pwm_list.append(pwm_row)
					if line[:4] == 'ID  ':
						id = line[4:-1]
				except:
					print 'Except: %s'%repr(sys.exc_info()[0])
					print line
					sys.exit(2)
			#output it if it's in profile_id_set
			if id in profile_id_set:
				profile_id_list.append(id)
				profile_id_outf.write('%s\t%s\n'%(len(profile_id_list), id))
				outf.write('>%s\n'%(len(pwm_list)))
				for pwm_row in pwm_list:
					outf.write('%.3f\t%.3f\t%.3f\t%.3f\n'%(pwm_row[0], pwm_row[1], pwm_row[2], pwm_row[3]))
			if self.report and block_no%500==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf, outf, profile_id_outf
		sys.stderr.write("Done\n")
	
	def run(self):
		profile_id_set = self.getProfileIdSet(self.profile_fname)
		self.transformMatrix(self.matrix_fname, self.output_fname, profile_id_set)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hp:m:o:ur", ["help", "debug", "report"])
	except:
		print __doc__
		sys.exit(2)
	
	profile_fname = None
	matrix_fname = None
	output_fname = None
	debug = 0
	report = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-p", ):
			profile_fname = arg
		elif opt in ("-m", ):
			matrix_fname = arg
		elif opt in ("-o", ):
			output_fname = arg
		elif opt in ("-u", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	if profile_fname and matrix_fname and output_fname:
		instance = TransformTRANSFACMatrixToCisMod(profile_fname, matrix_fname, output_fname, \
			debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
