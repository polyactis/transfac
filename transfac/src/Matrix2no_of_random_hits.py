#!/usr/bin/env python
"""
Usage: Matrix2no_of_random_hits.py -i input_dir -g GC_percentage [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, transfac(default)
	-i ...,	input_dir
	-g ...,	GC_percentage
	-o ...,	output_table, 'transfac.matrix2no_of_random_hits'(default)
	-c	commit the database transaction
	-b	enable debugging, no debug by default
	-r	report the progress(a number)
	-h, --help              show this help

Examples:
	Matrix2no_of_random_hits.py -i transfac/random_gc_95_out/ -c -g 0.95

Description:
	A program to fill in table matrix2no_of_random_hits. The input is TRANSFAC MATCH
		output, which is based on RandomSequenceGenerator.py's output.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt
from codense.common import db_connect

class Matrix2no_of_random_hits:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='transfac', 
		input_dir=None, GC_percentage=None, output_table='transfac.matrix2no_of_random_hits',\
		commit=0, debug=0, report=0):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_dir = input_dir
		self.GC_percentage = float(GC_percentage)
		self.output_table = output_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def submit2output_table(self, curs, output_table, seq_id, GC_percentage, mt_id2no_of_hits):
		for mt_id, no_of_hits in mt_id2no_of_hits.iteritems():
			curs.execute("insert into %s(mt_id, gc_perc, no_of_hits, seq_id)\
				values('%s', %s, %s, %s)"%(output_table, mt_id, GC_percentage, no_of_hits, seq_id))
	
	def parse_file(self, curs, input_fname, output_table, GC_percentage):
		sys.stderr.write("\t Parsing %s..."%os.path.basename(input_fname))
		inf = open(input_fname, 'r')
		mt_id2no_of_hits = {}
		seq_id = None
		for line in inf:
			if line[:10]=='Inspecting':
				if mt_id2no_of_hits and seq_id:	#1st time, both mt_id2no_of_hits and seq_id are nothing
					self.submit2output_table(curs, output_table, seq_id, GC_percentage, mt_id2no_of_hits)
				seq_id = int(line.split()[-1])	#\n is discarded automatically by split()
				mt_id2no_of_hits = {}	#clear the dictionary
			elif line[0] == ' ' and line[:6]!=' Total' and line[:6]!=' Frequ':	#the first character is blank, but exclude the end statistic part	
				ls = line[:-1].split('|')
				mt_id = ls[0].strip()	#remove spaces
				"""
				bs_disp_start_strand = ls[1].strip()
				bs_disp_start = int(bs_disp_start_strand[:-3])
				strand = bs_disp_start_strand[-2]
				core_similarity_score = ls[2]	#01-07-06
				matrix_similarity_score = ls[3]	#01-07-06
				sequence = ls[4].strip()
				#01-03-06
				bs_disp_end = bs_disp_start + len(sequence) - 1
				"""
				if mt_id not in mt_id2no_of_hits:
					mt_id2no_of_hits[mt_id] = 0
				mt_id2no_of_hits[mt_id] += 1
		
		if mt_id2no_of_hits and seq_id:	#don't forget the last sequence
			self.submit2output_table(curs, output_table, seq_id, GC_percentage, mt_id2no_of_hits)
		del inf
		sys.stderr.write("Done.\n")
	
	def run(self):
		files = os.listdir(self.input_dir)
		files.sort()
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		for input_fname in files:
			input_fname = os.path.join(self.input_dir, input_fname)
			self.parse_file(curs, input_fname, self.output_table, self.GC_percentage)
		
		if self.commit:
			curs.execute("end")
			
	
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:g:o:cbr", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'transfac'
	input_dir = None
	GC_percentage = None
	output_table = 'transfac.matrix2no_of_random_hits'
	commit = 0
	debug = 0
	report = 0

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-i",):
			input_dir = arg
		elif opt in ("-g",):
			GC_percentage = float(arg)
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
		
	if input_dir and output_table and GC_percentage:
		instance = Matrix2no_of_random_hits(hostname, dbname, schema, input_dir, GC_percentage, output_table, \
			commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
