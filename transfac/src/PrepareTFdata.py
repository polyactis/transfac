#!/usr/bin/env python
"""
Usage: PrepareTFdata.py [OPTIONS] -i INPUT_FILE -a GENE2ACC

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-i ...	INPUT_FILE
	-o ..., --output_table=...	output table. default is gene_id2tf
	-g ...	the two-letter abbreviation of organism
	-a ..., --acc_file=...	the gene2acc file
	-y ..., --type=...	input_file type, 1(refseq_id ...), 2(unigene_id refseq_id ...) 
	-u ..., --up_length=...	the upstream length, 1000(default)
	-m ..., --comment=...	the comment, ''(default)
	-n	output table needs to be created.
	-b, --debug	just test running the program, no daytime restriction
	-r, --report	report the progress (time before each query)
	-c, --commit	commit this database transaction
	-h, --help	show this help

Examples:
	PrepareTFdata.py -i /usr/local/research_data/tfdb/result.humanup2k.fptf
		-a  ncbi/gene/gene2accession  -u 2000 -m "5UTR, mask repeated as N" -r -c -g hs

Description:
	Parse Kangyu's TF result. Link its refseq_id to gene_id by looking
	up the gene2accession file.
	
"""

import psycopg, re, time, datetime, sys, urllib2, getopt, os, csv
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
sys.path += [os.path.join(os.path.expanduser('~/script/microarray/bin'))]
from MdbId2GeneId import MdbId2GeneId
from codense.common import db_connect, org2tax_id, org_short2long
from sets import Set

class PrepareTFdata:
	def __init__(self, hostname='zhoudb', dbname='mdb', input_filename=None, \
		output_table='gene_id2tf', organism=None, new_table=0, acc_file=None, \
		type=1, up_length=1000, comment='', debug=0, report=0, commit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.input_filename = input_filename
		self.output_table = output_table
		self.organism = organism
		self.new_table = int(new_table)
		self.acc_file = acc_file
		self.type = int(type)
		self.up_length = int(up_length)
		self.comment = comment
		self.debug = int(debug)
		self.report = int(report)
		self.commit = int(commit)
	
	def parse_input_filename(self, curs, input_filename, output_table, acc2gene_id, tax_id, up_length, comment, organism, type=1):
		"""
		09-07-05
		09-08-05 fix a bug and improve the status reporting.
		"""
		sys.stderr.write("Parsing %s ...\n"%input_filename)
		inf = open(input_filename, 'r')
		counter = 0
		no_of_mapping_fails = 0
		for line in inf:
			ls = line.split()
			if type==2:
				ls.pop(0)	#throw away the unigene_id
			refseq_id = ls.pop(0)
			if refseq_id[-1] == '*':
				refseq_id = refseq_id[:-1]	#throw away '*'
			ls.pop(0)	#throw away the '|'
			
			key = (refseq_id.upper(), tax_id)
			if key in acc2gene_id:
				while len(ls)>0:
					tf = ls.pop(0)
					freq = int(ls.pop(0))
					for gene_id in acc2gene_id[key]:
						row = [gene_id, tf, freq, up_length, organism, comment]
						self.submit2table(curs, output_table, row)
			else:
				no_of_mapping_fails+=1
			counter+=1
			if self.report and counter%1000==0:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		if self.report:
			sys.stderr.write('%s%s'%('\x08'*20, counter))
		del inf
		sys.stderr.write("(%s refseq_id no gene_id mapped) Done.\n"%no_of_mapping_fails)
	
	def submit2table(self, curs, table_name, row):
		"""
		09-08-05
		"""
		if row[5] =='':
			curs.execute("insert into graph.%s(gene_id, tf, freq, up_length, organism)\
				values('%s', '%s', %s, %s, '%s')"%\
				(table_name, row[0], row[1], row[2], row[3], row[4]))
		else:
			curs.execute("insert into graph.%s(gene_id, tf, freq, up_length, organism,\
				comment) values('%s', '%s', %s, %s, '%s', '%s')"%\
				(table_name, row[0], row[1], row[2], row[3], row[4], row[5]))
	
	def create_output_table(self, curs, output_table):
		"""
		09-08-05
		"""
		sys.stderr.write("Creating table %s ..."%output_table)
		curs.execute("create table graph.%s(\
				gene_id	varchar,\
				tf	varchar,\
				freq	integer,\
				up_length	integer,\
				organism	varchar,\
				comment	varchar)"%output_table)
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		09-08-05
			
			--db_connect()
			--org_short2long()
			--org2tax_id()
			--setup_acc2gene_id()
			if self.new_table
				--create_output_table()
			--parse_input_filename()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname)
		long_organism = org_short2long(self.organism)
		tax_id_set = Set([org2tax_id(long_organism)])
		
		MdbId2GeneId_instance = MdbId2GeneId()
		acc2gene_id = MdbId2GeneId_instance.setup_acc2gene_id(self.acc_file, tax_id_set)
		if self.new_table:
			self.create_output_table(curs, self.output_table)
		self.parse_input_filename(curs, self.input_filename, self.output_table,  acc2gene_id,\
			org2tax_id(long_organism), self.up_length, self.comment, long_organism, self.type)
		if self.commit:
			curs.execute("end")
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "output_table=", "acc_file=", "type= ", \
		"up_length=", "comment=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:i:o:g:a:y:u:m:nbrch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	input_filename = None
	output_table = 'gene_id2tf'
	organism = None
	acc_file = None
	type = 1
	up_length = 1000
	comment = ''
	new_table = 0
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
		elif opt in ("-i"):
			input_filename = arg
		elif opt in ("-o", "--output_table"):
			output_table = arg
		elif opt in ("-g"):
			organism = arg
		elif opt in ("-a", "--acc_file"):
			acc_file = arg
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("-u", "--up_length"):
			up_length = int(arg)
		elif opt in ("-m", "--comment"):
			comment = arg
		elif opt in ("-n"):
			new_table = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1

	if input_filename and acc_file and organism:
		instance = PrepareTFdata(hostname, dbname, input_filename, output_table, \
			organism, new_table, acc_file, type, up_length, comment, debug, report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
