#!/usr/bin/env python
"""
Usage: LinkFactorID2EntrezGeneID.py -a

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, transfac(default)
	-f ...,	factor_table, transfac.factor(default)
	-a ...,	the gene2acc file
	-o ...,	output table, transfac.tf_acc2entrezgene_id(default)
	-p ...,	another output table to store raw linking results, transfac.tf_acc2entrezgene_id_raw(default)
	-c,	commit the database transaction
	-b,	enable debug flag
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	~/script/annot/bin/LinkFactorID2EntrezGeneID.py -a /usr/local/research_data/ncbi/gene/gene2accession

Description:
	Program to link tf_acc of TRANSFAC's factors to Entrez Gene ID.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv, re
from codense.common import db_connect, get_tax_id_from_org
from sets import Set

class LinkFactorID2EntrezGeneID:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='transfac', factor_table='transfac.factor', \
		gene2acc_file=None, output_table='transfac.tf_acc2entrezgene_id', \
		raw_output_table='transfac.tf_acc2entrezgene_id_raw', commit=0, debug=0, report=0):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.factor_table = factor_table
		self.gene2acc_file = gene2acc_file
		self.output_table = output_table
		self.raw_output_table = raw_output_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
		
	
	def get_factor_info(self, curs, factor_table):
		"""
		01-31-06
			one typical external_database_links looks like:
			{"TRANSPATH: MO000024633","EMBL: X12549; DMASCT3","SWISSPROT: P09774; AST3_DROME","PIR: S01165; S01165","FLYBASE: FBgn0002561; l(1)sc"}
			
			only focus on EMBL or SWISSPROT 
		"""
		sys.stderr.write("Getting factor info...")
		curs.execute("select tf_acc, organism, external_database_links from %s \
			where external_database_links is not null"%factor_table)
		rows = curs.fetchall()
		acc_tax_id2tf_acc = {}
		organism2tax_id = {}
		for row in rows:
			tf_acc, organism, external_database_links = row
			if organism in organism2tax_id:
				tax_id = organism2tax_id[organism]
			else:
				tax_id = get_tax_id_from_org(curs, organism)
				organism2tax_id[organism] = tax_id
			if tax_id:	#not some weird name that don't have tax_id
				external_database_links = external_database_links[2:-2].split('","')
				for external_database_link in external_database_links:
					xdb_name_acc_ls = external_database_link.split(':')
					xdb_name = xdb_name_acc_ls[0]
					if xdb_name=='EMBL' or xdb_name=='SWISSPROT':
						xdb_acc = xdb_name_acc_ls[1]
						xdb_acc = xdb_acc.split(';')[0].strip()
						key = (xdb_acc.upper(), tax_id)
						acc_tax_id2tf_acc[key] = tf_acc
		sys.stderr.write("Done.\n")
		return acc_tax_id2tf_acc
	
	def setup_acc2gene_id(self, gene2acc_file, acc_tax_id2tf_acc):
		"""
		01-31-06
			copied from MdbId2GeneId.py
		"""
		sys.stderr.write("Setting up acc2gene_id...")
		p_gb_acc_version = re.compile(r'(^\w+)\.\d+')
		tf_acc2gene_id_bridge_acc_ls = {}
		reader = csv.reader(open(gene2acc_file, 'r'), delimiter ='\t')
		for row in reader:
			tax_id, gene_id, status, rna_nuc_acc_ver, rna_nuc_gi, prot_acc_ver,\
			prot_gi, genomic_nuc_acc_ver, genomic_nuc_gi, start, end, orientation = row
			tax_id = int(tax_id)	#integer
			gene_id = int(gene_id)	#integer
			acc_ver_to_check = [rna_nuc_acc_ver, rna_nuc_gi, prot_acc_ver,\
				prot_gi, genomic_nuc_acc_ver, genomic_nuc_gi]
			for acc_ver in acc_ver_to_check:
				if acc_ver!='-':
					key = (acc_ver.upper(), tax_id)
					if key in acc_tax_id2tf_acc:
						tf_acc=  acc_tax_id2tf_acc[key]
						if tf_acc not in tf_acc2gene_id_bridge_acc_ls:
							tf_acc2gene_id_bridge_acc_ls[tf_acc] = []
						tf_acc2gene_id_bridge_acc_ls[tf_acc].append((gene_id, key[0]))
					p_gb_acc_version_result = p_gb_acc_version.search(acc_ver)
					if p_gb_acc_version_result:
						accession = p_gb_acc_version_result.groups()[0]
						key = (accession.upper(), tax_id)
						if key in acc_tax_id2tf_acc:
							tf_acc=  acc_tax_id2tf_acc[key]
							if tf_acc not in tf_acc2gene_id_bridge_acc_ls:
								tf_acc2gene_id_bridge_acc_ls[tf_acc] = []
							tf_acc2gene_id_bridge_acc_ls[tf_acc].append((gene_id, key[0]))
		del reader
		sys.stderr.write("Done\n")
		return tf_acc2gene_id_bridge_acc_ls
	
	def submit_raw_result(self, curs, tf_acc2gene_id_bridge_acc_ls, raw_output_table):
		"""
		01-31-06
			Watch: tf_accs with >1 gene_id linked are also submitted.
			Some genomic_nuc_acc_ver leads to numerous gene_id linked to the tf_acc,
			which is apparently faulty.
		"""
		sys.stderr.write("Submitting to %s..."%raw_output_table)
		tf_acc2entrezgene_id_set = {}
		for tf_acc, gene_id_bridge_acc_ls in tf_acc2gene_id_bridge_acc_ls.iteritems():
			tf_acc2entrezgene_id_set[tf_acc] = Set()
			for gene_id_bridge_acc in gene_id_bridge_acc_ls:
				curs.execute("insert into %s(tf_acc, gene_id, bridge_acc) values('%s', %s, '%s')"%\
					(raw_output_table, tf_acc, gene_id_bridge_acc[0], gene_id_bridge_acc[1]))
				tf_acc2entrezgene_id_set[tf_acc].add(gene_id_bridge_acc[0])
		sys.stderr.write("Done\n")
		return tf_acc2entrezgene_id_set
	
	def submit_result(self, curs, tf_acc2entrezgene_id_set, output_table):
		"""
		01-31-06
			disregard those tf_acc's with >1 gene_id
		"""
		sys.stderr.write("Submitting to %s..."%output_table)
		for tf_acc, entrezgene_id_set in tf_acc2entrezgene_id_set.iteritems():
			if len(entrezgene_id_set)==1:	#only those with 1 gene_id
				curs.execute("insert into %s(tf_acc, gene_id) values('%s', %s)"%(output_table, tf_acc, entrezgene_id_set.pop()))
		sys.stderr.write("Done\n")
	
	def run(self):
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		acc_tax_id2tf_acc = self.get_factor_info(curs, self.factor_table)
		tf_acc2gene_id_bridge_acc_ls = self.setup_acc2gene_id(self.gene2acc_file, acc_tax_id2tf_acc)
		tf_acc2entrezgene_id_set = self.submit_raw_result(curs, tf_acc2gene_id_bridge_acc_ls, self.raw_output_table)
		self.submit_result(curs, tf_acc2entrezgene_id_set, self.output_table)
		
		if self.commit:
			curs.execute("end")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:f:a:o:p:cbr", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'transfac'
	factor_table = 'transfac.factor'
	gene2acc_file = None
	output_table = 'transfac.tf_acc2entrezgene_id'
	raw_output_table = 'transfac.tf_acc2entrezgene_id_raw'
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
		elif opt in ("-f",):
			factor_table = arg
		elif opt in ("-a",):
			gene2acc_file = arg
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-p",):
			raw_output_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
		
	if factor_table and gene2acc_file and output_table and raw_output_table:
		instance = LinkFactorID2EntrezGeneID(hostname, dbname, schema, factor_table, \
			gene2acc_file, output_table, raw_output_table, commit,  debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
