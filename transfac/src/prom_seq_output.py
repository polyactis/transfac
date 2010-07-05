#!/usr/bin/env python
"""
Usage: prom_seq_output.py [OPTIONS] -f FOLDER

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, transfac(default), or sequence
	-f ...	the folder to store the sequence files
	-p ...	output file prefix, numbered from 0, prom_seq_output(default)
	-s ...	number of sequences per file, 10,000(default)
	-g ...,	organism, hs(default)
	-y ...,	running type, 1(default, from table prom_seq), 2(from table raw_sequence, set schema accordingly)
	-b, --debug	just test running the program, no daytime restriction
	-r, --report	report the progress (time before each query)
	-c, --commit	commit this database transaction(IGNORE)
	-h, --help	show this help

Examples:
	prom_seq_output.py -f flanking_seq
	prom_seq_output.py -k sequence -y2 -f genome_seq

Description:
	Program to output sequences from prom_seq in chunks.
	
	or from table raw_sequence
	
"""

import psycopg, sys, getopt, os
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
sys.path += [os.path.join(os.path.expanduser('~/script/microarray/bin'))]
from codense.common import db_connect, org_short2long, org2tax_id

class prom_seq_output:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='transfac', folder=None, \
		prefix='prom_seq_output', size=10000, organism='hs', running_type=1, debug=0, report=0, commit=0):
		"""
		2006-11-27
			add running_type
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.folder = folder
		self.prefix = prefix
		self.size = int(size)
		self.organism = org_short2long(organism)
		self.running_type = int(running_type)
		self.debug = int(debug)
		self.report = int(report)
		self.commit = int(commit)
	
	def run(self):
		"""
		11-15-05
			correct a bug related to self.size
		2006-08-27
			if sequence is empty, ignore it.
		2006-11-27
			add running_type
		"""
		if not os.path.isdir(self.folder):
			os.makedirs(self.folder)
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		if self.running_type==1:
			curs.execute("DECLARE crs CURSOR FOR SELECT id, sequence from prom_seq \
				where sequence is not null and strpos(chromosome, 'random')=0 and organism='%s'"%self.organism)
				#09-14-05	not null sequence and no 'random' in chromosome
		elif self.running_type==2:
			curs.execute("DECLARE crs CURSOR FOR SELECT r.id, r.sequence from sequence.raw_sequence r, sequence.annot_assembly a\
				where r.acc_ver=a.acc_ver and a.tax_id=%s"%org2tax_id(self.organism))
				#2006-11-27 from a specific tax_id, the later condition will guarantee that sequence is not empty
		else:
			sys.stderr.write("Unsupported running_type: %s\n"%self.running_type)
			sys.exit(3)
		curs.execute("fetch %s from crs"%self.size)
		rows = curs.fetchall()
		counter = 0
		sys.stderr.write("Starting to output...\n")
		while rows:
			output_file = os.path.join(self.folder, '%s%s'%(self.prefix, counter))
			of = open(output_file, 'w')
			for row in rows:
				id, sequence = row
				if sequence:
					of.write('>%s\n%s\n'%(id,sequence))
			del of
			counter += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			curs.execute("fetch %s from crs"%self.size)
			rows = curs.fetchall()
		del conn, curs
		sys.stderr.write("Done.\n")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:f:p:s:g:y:brch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'transfac'
	folder = None
	prefix = 'prom_seq_output'
	size = 10000
	organism = 'hs'
	running_type = 1
	debug = 0
	report = 0
	commit = 0
	
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
		elif opt in ("-f",):
			folder = arg
		elif opt in ("-p",):
			prefix = arg
		elif opt in ("-s",):
			size = int(arg)
		elif opt in ("-g",):
			organism = arg
		elif opt in ("-y",):
			running_type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1

	if folder and hostname and dbname and schema and prefix and size:
		instance = prom_seq_output(hostname, dbname, schema, folder, \
			prefix, size, organism, running_type, debug, report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
