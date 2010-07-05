#!/usr/bin/env python
"""
Usage: TRANSFACHits2gene_id2mt_no.py -i input_fname -g organism [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-i ...,	input_fname
	-m ...,	the matrix table(transfac.matrix, default)
	-q ...,	the prom_seq_table(transfac.prom_seq, default)
	-o ...,	the output table(graph.gene_id2mt_no_p_value, default)
	-p ...,	p_value_cut_off, 0.001(default)
	-g ...,	organism
	-c	commit
	-b	debug version.
	-r	enable report flag
	-h, --help	Display the usage infomation.
	
Examples:
	TRANSFACHits2gene_id2mt_no.py -i transfac/prom_seq_hs_match_out.hits.pvalue.data
		-g hs -c -r
	
Description:
	Construct gene_id2mt_no_p_value from AnalyzeTRANSFACHits.py's output.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv
from codense.common import db_connect, org_short2long, get_mt_id2no, get_tax_id_from_org
from sets import Set

class TRANSFACHits2gene_id2mt_no:
	def __init__(self, hostname='zhoudb', dbname='graphdb', input_fname=None,\
		matrix_table='transfac.matrix', prom_seq_table='transfac.prom_seq', \
		output_table='graph.gene_id2mt_no_p_value', p_value_cut_off=0.001, organism='hs', \
		commit=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.input_fname = input_fname
		self.matrix_table = matrix_table
		self.prom_seq_table = prom_seq_table
		self.output_table = output_table
		self.p_value_cut_off = float(p_value_cut_off)
		self.organism = org_short2long(organism)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_prom_id2gene_id(self, curs, prom_seq_table, organism):
		sys.stderr.write("Getting prom_id2gene_id...\n")
		curs.execute("DECLARE crs_p2g CURSOR FOR select id, prom_acc from %s where organism='%s'"%(prom_seq_table, organism))
		curs.execute("fetch 10000 from crs_p2g")
		rows = curs.fetchall()
		prom_id2gene_id = {}
		counter = 0
		while rows:
			for row in rows:
				id, prom_acc = row
				prom_id2gene_id[id] = int(prom_acc)
				counter += 1
			if self.report:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
			curs.execute("fetch 10000 from crs_p2g")
			rows = curs.fetchall()
		if self.report:
			sys.stderr.write("%s%s"%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		return prom_id2gene_id
	
	def submit_hits(self, curs, row, tax_id, output_table):
		gene_id, sequence_length, mt_no, gc_perc, no_of_hits, p_value = row
		curs.execute("INSERT INTO %s(gene_id, mt_no, seq_length, no_of_hits, gc_perc, p_value, tax_id)\
			values(%s, %s, %s, %s, %s, %s, %s)"%(output_table, gene_id, mt_no, sequence_length, \
			no_of_hits, gc_perc, p_value, tax_id))
	
	
	def parse_input_fname(self, curs, input_fname, p_value_cut_off, prom_id2gene_id, mt_id2no, tax_id, output_table):
		sys.stderr.write("Parsing %s...\n"%os.path.basename(input_fname))
		reader = csv.reader(open(input_fname,'r'), delimiter='\t')
		counter = 0
		real_counter = 0
		for row in reader:
			seq_id, sequence_length, mt_id, gc_perc, no_of_hits, p_value = row
			p_value = float(p_value)
			if p_value<=p_value_cut_off:
				seq_id = int(seq_id)
				gene_id = prom_id2gene_id[seq_id]
				mt_no = mt_id2no[mt_id]
				new_row = [gene_id, sequence_length, mt_no, gc_perc, no_of_hits, p_value]
				self.submit_hits(curs, new_row, tax_id, output_table)
				real_counter += 1
			counter += 1
			if self.report and counter%10000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*20, counter, real_counter))
		if self.report and counter%10000==0:
			sys.stderr.write("%s%s\t%s"%('\x08'*20, counter, real_counter))
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		02-01-06
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname)
		tax_id = get_tax_id_from_org(curs, self.organism)
		mt_id2no = get_mt_id2no(curs, self.matrix_table)
		prom_id2gene_id = self.get_prom_id2gene_id(curs, self.prom_seq_table, self.organism)
		
		self.parse_input_fname(curs, self.input_fname, self.p_value_cut_off, prom_id2gene_id, mt_id2no, tax_id, self.output_table)
		if self.commit:
			curs.execute("end")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:i:m:q:o:p:g:brch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	input_fname = ''
	matrix_table = 'transfac.matrix'
	prom_seq_table= 'transfac.prom_seq'
	output_table = 'graph.gene_id2mt_no_p_value'
	p_value_cut_off = 0.001
	organism = ''
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
		elif opt in ("-i", ):
			input_fname = arg
		elif opt in ("-m",):
			matrix_table = arg
		elif opt in ("-q",):
			prom_seq_table = arg
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-p",):
			p_value_cut_off = float(arg)
		elif opt in ("-g",):
			organism = arg
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
		elif opt in ("-c",):
			commit = 1

	if input_fname and organism:
		instance = TRANSFACHits2gene_id2mt_no(hostname, dbname, input_fname, matrix_table, \
			prom_seq_table, output_table, p_value_cut_off, organism, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
