#!/usr/bin/env python
"""
Usage: binding_site2gene_id2mt_no.py-a GENE2ACCESSION [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-i ...,	the binding_site table(transfac.binding_site, default)
	-m ...,	the matrix table(transfac.matrix, default)
	-p ...,	the prom_seq_table(transfac.prom_seq, default)
	-o ...,	the output table(graph.gene_id2mt_no, default)
	-t ...,	top number of genes selected , 2000(default)
	-a ...,	the prom_acc to gene_id linking file, gene2accession or gene_info
	-g ...,	organism, hs(default)
	-c	commit
	-b	debug version.
	-r	enable report flag
	-h, --help	Display the usage infomation.
	
Examples:
	binding_site2gene_id2mt_no.py -a /usr/local/research_data/ncbi/gene/gene2accession -c -r
	binding_site2gene_id2mt_no.py -a /usr/local/research_data/ncbi/gene/gene2accession -g mm -c -r
	# a total different version, try on the at_tfdb
	binding_site2gene_id2mt_no.py -a /usr/local/research_data/ncbi/gene/gene_info
		-i at_tfdb.binding_site -m at_tfdb.matrix  -p at_tfdb.prom_seq -t 763508 -g at -c -r
Description:
	Construct gene_id2mt_no from binding_site. top_number is specific for one organism.
	So one organism one run.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/microarray/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/microarray/bin')))
import sys, os, getopt, csv
from MdbId2GeneId import MdbId2GeneId
from codense.common import db_connect, org2tax_id, org_short2long, get_mt_id2no
from heapq import heappush, heappop
from sets import Set
from threading import Thread

class get_gene_id_dict(Thread):
	"""
	09-30-05
		a thread to speed up the main program
	"""
	def __init__(self, input_file, tax_id_set):
		Thread.__init__(self)
		self.input_file = input_file
		self.tax_id_set = tax_id_set
		self.gene_id_dict = None
		
	def run(self):
		MdbId2GeneId_instance = MdbId2GeneId()
		filename = os.path.basename(self.input_file)
		if filename=='gene2accession':
			self.gene_id_dict = MdbId2GeneId_instance.setup_acc2gene_id(self.input_file, self.tax_id_set)
		elif filename=='gene_info':
			self.gene_id_dict = MdbId2GeneId_instance.setup_name2gene_id(self.input_file, self.tax_id_set)
		else:
			sys.stderr.write("%s is neither gene2accession nor gene_info. Abort.\n"%self.acc_file)
			sys.exit(3)
		del MdbId2GeneId_instance

class attr_of_mt_no:
	"""
	09-18-05
		data structure to store [matrix_similarity_score, prom_acc, id]
		The idea is to store them while keep the number of distinct prom_acc's under top_number
	"""
	def __init__(self, top_number, debug=0):
		self.top_number = top_number
		self.debug = int(debug)
		self.hq_list = []
		self.acc2id_set = {}
	
	def remove_lowest_row(self):
		"""
		09-18-05
			remove the row with lowest_score
		"""
		lowest_score, old_prom_acc, old_id = heappop(self.hq_list)		
		self.acc2id_set[old_prom_acc].remove(old_id)
		if len(self.acc2id_set[old_prom_acc]) == 0:	#remove the old_prom_acc if there's no more old_prom_acc in the hq_list
			del self.acc2id_set[old_prom_acc]
	
	def add_new_row(self, matrix_similarity_score, prom_acc, id):
		heappush(self.hq_list, [matrix_similarity_score, prom_acc, id])
		if prom_acc not in self.acc2id_set:	#add the new prom_acc
			self.acc2id_set[prom_acc] = Set()
		self.acc2id_set[prom_acc].add(id)
	
	def consume_new_row(self, matrix_similarity_score, prom_acc, id):
		if len(self.acc2id_set)==self.top_number:	#no. of distinct accs exceeds the limit
			lowest_score = self.hq_list[0][0]
			if lowest_score<matrix_similarity_score:	#if only the current score is higher than the lowest_score
				self.remove_lowest_row()
				self.add_new_row(matrix_similarity_score, prom_acc, id)
				while len(self.acc2id_set)>self.top_number:	#IF the old_prom_acc has a higher score match AND the new prom_acc is distinct
					self.remove_lowest_row()	#remove the row with lowest_score until no. of distinct accs is under top_number.
					
		elif len(self.acc2id_set)<self.top_number:
			self.add_new_row(matrix_similarity_score, prom_acc, id)
		else:
			sys.stderr.write('%s\t%s\t%s\n'%(matrix_similarity_score, prom_acc, id))
			sys.stderr.write("Heap Queue's size exceeds the limit. Something impossible. Exit\n")
			sys.exit(2)
	
	def release(self):
		"""
		09-16-05
			throw away hq_list to release memory
		"""
		self.lowest_score = self.hq_list[0][0]
		del self.hq_list
		
class binding_site2gene_id2mt_no:
	def __init__(self, hostname='zhoudb', dbname='graphdb', input_table='transfac.binding_site',\
		matrix_table='transfac.matrix', prom_seq_table='transfac.prom_seq', \
		output_table='graph.gene_id2mt_no', top_number=2000, acc_file=None, organism='hs', \
		commit=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.input_table = input_table
		self.matrix_table = matrix_table
		self.prom_seq_table = prom_seq_table
		self.output_table = output_table
		self.top_number = int(top_number)
		self.acc_file = acc_file
		self.organism = org_short2long(organism)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_mt_no2matches(self, curs, input_table, prom_seq_table, top_number, organism, mt_id2no):
		"""
		09-18-05
			use heapq to keep each mt_no's associated matches <= top_number
		09-30-05
			make prom_seq_table explicit
		"""
		sys.stderr.write("Getting mt_no2matches from %s...\n"%input_table)
		mt_no2matches = {}
		curs.execute("DECLARE crs CURSOR FOR select b.mt_id, b.matrix_similarity_score, p.prom_acc, \
			b.id, p.organism from %s b, %s p where p.id=b.prom_id"%(input_table, prom_seq_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				mt_id, matrix_similarity_score, prom_acc, id, org	= row
				if org!=organism:
					continue
				mt_no = mt_id2no[mt_id]
				if mt_no not in mt_no2matches:
					mt_no2matches[mt_no] = attr_of_mt_no(top_number, self.debug)
				mt_no2matches[mt_no].consume_new_row(matrix_similarity_score, prom_acc, id)
				counter +=1
				
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			if self.debug and counter == 100000:
				break
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		
		for mt_no, unit in mt_no2matches.iteritems():	#finally, release some memory
			unit.release()
		sys.stderr.write("Done.\n")
		return mt_no2matches
	
	def submit2output_table(self, curs, output_table, mt_no, gene_id2id_list, tax_id):
		"""
		09-18-05
		"""
		for gene_id, id_list in gene_id2id_list.iteritems():
			id_list.sort()
			id_list = '{' + repr(id_list)[1:-1] + '}'
			curs.execute("insert into %s(gene_id, mt_no, binding_site_id_list, tax_id)\
				values('%s', '%s', '%s', %s)"%(output_table, gene_id, mt_no, id_list, tax_id))
		
	def dump2output_table(self, curs, output_table, mt_no2matches, acc2gene_id, tax_id):
		"""
		09-18-05
			one gene_id could have >=2 acc's, so need to transfer data to gene_id2id_list
		09-29-05
			check no-version-form of the acc to link to gene_id
		"""
		sys.stderr.write("Dumping to %s...\n"%output_table)
		counter = 0
		perfect_matrices = []
		MdbId2GeneId_instance = MdbId2GeneId()
		for mt_no, unit in mt_no2matches.iteritems():
			if unit.lowest_score==1.0:
				perfect_matrices.append(mt_no)
				continue
			gene_id2id_list = {}
			for acc, id_set in unit.acc2id_set.iteritems():
				key = (acc.upper(), tax_id)
				gene_id_list = acc2gene_id.get(key)
				if not gene_id_list:
					p_gb_acc_version_result = MdbId2GeneId_instance.p_gb_acc_version.search(acc)
					if p_gb_acc_version_result:
						acc_no_ver = p_gb_acc_version_result.groups()[0]
						key = (acc_no_ver.upper(), tax_id)
						gene_id_list = acc2gene_id.get(key)
						if not gene_id_list:	#still nothing
							if self.report:
								sys.stderr.write("%s link-failure.\n"%acc)
							continue
					else:	#no version
						if self.report:
							sys.stderr.write("%s link-failure.\n"%acc)
						continue
				for gene_id in gene_id_list:
					if gene_id not in gene_id2id_list:
						gene_id2id_list[gene_id] = []
					for id in id_set:
						gene_id2id_list[gene_id].append(id)
						counter +=1
			self.submit2output_table(curs, output_table, mt_no, gene_id2id_list, tax_id)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		sys.stderr.write("Perfect_matrices skipped: %s\n"%repr(perfect_matrices))
	
	def run(self):
		"""
		09-18-05
		09-30-05
			way of calling get_mt_id2no() changed
			get gene_id_dict is a thread now, to speed up
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname)
		tax_id = org2tax_id(self.organism)
		tax_id_set = Set([tax_id])
		get_gene_id_dict_instance = get_gene_id_dict(self.acc_file, tax_id_set)
		get_gene_id_dict_instance.start()
		mt_id2no = get_mt_id2no(curs, self.matrix_table)
		mt_no2matches = self.get_mt_no2matches(curs, self.input_table, self.prom_seq_table, self.top_number, self.organism, mt_id2no)
		get_gene_id_dict_instance.join()	#must wait it to finish before going on, need gene_id_dict
		self.dump2output_table(curs, self.output_table, mt_no2matches, get_gene_id_dict_instance.gene_id_dict, tax_id)
		if self.commit:
			curs.execute("end")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:i:m:p:o:t:a:g:brch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	input_table = 'transfac.binding_site'
	matrix_table = 'transfac.matrix'
	prom_seq_table= 'transfac.prom_seq'
	output_table = 'graph.gene_id2mt_no'
	top_number = 2000
	acc_file = None
	organism = 'hs'
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
			input_table = arg
		elif opt in ("-m"):
			matrix_table = arg
		elif opt in ("-p"):
			prom_seq_table = arg
		elif opt in ("-o"):
			output_table = arg
		elif opt in ("-t"):
			top_number = int(arg)
		elif opt in ("-a"):
			acc_file = arg
		elif opt in ("-g"):
			organism = arg
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
		elif opt in ("-c"):
			commit = 1

	if acc_file:
		instance = binding_site2gene_id2mt_no(hostname, dbname, input_table, matrix_table, \
			prom_seq_table, output_table, top_number, acc_file, organism, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
