#!/usr/bin/env python
"""
Usage: transfacdb.py [OPTIONS] -i INPUT_FILE -y PARSER_TYPE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, transfac(default)
	-i ...	INPUT_FILE
	-p ...	sequence_type(1 upstream(default), 2 5UTR, 3 3UTR, 4 downstream)
		only for type=1
	-g ...	the two-letter abbreviation of organism, hs(default)
		only for type =1
	-y ...	parser type(1, default)
	-b, --debug	just test running the program, no daytime restriction
	-r, --report	report the progress (time before each query)
	-c, --commit	commit this database transaction
	-h, --help	show this help

Examples:
	./transfacdb.py -i rn_up2k.seq -y 1 -p1 -c -r -g rn
	./transfacdb.py -i ../data/factor.dat -y2  -c -r
	./transfacdb.py -i ../data/matrix.dat -y3 -c -r

Description:
	program to parse several kinds of files and dump the data into schema transfac
	1. prom_seq	2. factor	3.matrix	4.binding_site
	5.reference	6.site	7.cell	8.fragment	9.gene
	10.class	11. binding_site_easy_parse(output of TopMatchResult.py)
	
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
import psycopg, sys, getopt, csv, cStringIO, re
from codense.common import db_connect, org_short2long
from MdbId2UnigeneId import unigene_data_block_iterator

class fasta_block_iterator:
	'''
	09-13-05
		fasta format iterator
		a little bit tricky, '>', the block starter is used as a tokenizer
	2006-09-01
		it seems 'for line in self.inf' doesn't work on hpc-cmb.
		check https://dl403k-1.cmb.usc.edu/log/hpc-cmb
	'''
	def __init__(self, inf):
		self.inf = inf
		self.block = ''
		self.previous_line = ''
	def __iter__(self):
		return self
	def next(self):
		self.read()
		return self.block
	def read(self):
		self.block = self.previous_line	#don't forget the starting line
		line = self.inf.readline()
		while(line):
			if line[0]=='>':
				self.previous_line = line
				if self.block:	#not the first time
					break
				else:	#first time to read the file, block is still empty
					self.block += line
			else:
				self.block += line
			line = self.inf.readline()
		if self.block==self.previous_line:
			raise StopIteration

class match_block_iterator:
	'''
	09-14-05
		iterator for match's output file
	2006-09-01
		same modification as fasta_block_iterator
	'''
	def __init__(self, inf):
		self.inf = inf
		self.block = ''
		self.previous_line = ''
	def __iter__(self):
		return self
	def next(self):
		self.read()
		return self.block
	def read(self):
		self.block = self.previous_line	#don't forget the starting line
		line = self.inf.readline()
		while(line):
			if line[:10]=='Inspecting':
				self.previous_line = line
				if self.block:	#not the first time
					break
				else:	#first time to read the file, block is still empty
					self.block += line
			elif line[0] == ' ' and line[:6]!=' Total' and line[:6]!=' Frequ':	#the first character is blank, but exclude the end statistic part
				self.block += line
			line = self.inf.readline()
		if self.block==self.previous_line:	#no stuff in the file anymore, the previous for loop is skipped
			raise StopIteration

class transfac_embl_block:
	"""
	09-15-05
		a data structure to hold all information of a transfac_embl_block
	"""
	def __init__(self):
		self.acc = ''
		self.basis = ''
		self.binding_region = '{}'
		self.cell_description = ''
		self.cell_factor_source = ''
		self.cell_negative = ''
		self.cell_positive = ''
		self.cell_source = '{}'
		self.chromosome = ''
		self.class_acc = ''
		self.class_tree_id = ''
		self.class_struc = ''
		self.comment = ''
		self.consensus = ''
		self.date = ''
		self.description = ''
		self.element = ''
		self.expr_pattern = '{}'
		self.external_database_links = '{}'
		self.factors = '{}'
		self.feature = '{}'
		self.func_feature = '{}'
		self.gene_acc = ''
		self.gene_region = ''
		self.homolog = ''
		self.id = ''
		self.interacting_partners = '{}'
		self.matrices = '{}'
		self.method = '{}'
		self.name = ''
		self.organism = ''
		self.position_start = 2000000
		self.position_end = 2000000
		self.position_ref_point = ''
		self.prom_classification = ''
		self.pwm = '{}'
		self.reference = '{}'
		self.ref_external_link = ''
		self.ref_authors = ''
		self.ref_title = ''
		self.ref_journal = ''
		self.regulation = ''
		self.sequence = ''
		self.sequence_source = ''
		self.short_description = ''
		self.site_in_matrix_accs = '{}'
		self.site_in_matrix_list = []
		self.sites = '{}'
		self.size = ''
		self.struc_feature = '{}'
		self.synonyms = ''
		self.tf_name = ''
		self.type = ''


class transfacdb:
	"""
	09-13-05
	"""
	def __init__(self, hostname='zhoudb', dbname='mdb', schema='transfac', inputfile=None, \
		output_table=None, sequence_type=1, organism='hs', type=1, debug=0, report=0, commit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputfile = inputfile
		self.output_table = output_table
		self.sequence_type = int(sequence_type)
		self.organism = org_short2long(organism)
		self.type = int(type)
		self.debug = int(debug)
		self.report = int(report)
		self.commit = int(commit)
		self.parser_dict = {1: self.prom_seq_parse,
			2: self.factor_parse,
			3: self.matrix_parse,
			4: self.binding_site_parse,
			5: self.reference_parse,
			6: self.site_parse,
			7: self.cell_parse,
			8: self.fragment_parse,
			9: self.gene_parse,
			10: self.class_parse,
			11: self.binding_site_easy_parse}
		self.output_table_dict = {1: 'prom_seq',
			2: 'factor',
			3: 'matrix',
			4: 'binding_site',
			5: 'reference',
			6: 'site',
			7: 'cell',
			8: 'fragment',
			9: 'gene',
			10: 'class',
			11: 'binding_site'}
		self.pwm_line_pattern = re.compile(r'\d\d  ')
		
	def submit2prom_seq(self, curs, table_name, row):
		curs.execute("insert into %s(prom_acc, chromosome, strand, prom_genome_start, prom_genome_end, organism, prom_type_id, sequence)\
			values('%s', '%s', '%s', %s, %s, '%s', %s, '%s')"%\
			(table_name, row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7]))
	
	def parse_fasta_first_line(self, header_line):
		"""
		header_line:
		>rn3_knownGene_AY389467 range=chr1:1304646-1306645 5'pad=0 3'pad=0 revComp=FALSE strand=+ repeatMasking=lower
		or
		>hg17_knownGene_NM_021983 range=chr6_hla_hap2:18514-18584 5'pad=0 3'pad=0 revComp=TRUE strand=- repeatMasking=lower
		"""
		ls = header_line[1:-1].split(' ')
		acc = ls[0]
		acc = acc.split('_')[2:]	#-> ['NM', '021983']
		acc = '_'.join(acc)	#-> 'NM_021983'
		range =ls[1]
		range = range.split('=')	#-> ['range', 'chr1:1304646-1306645']
		range = range[1].split(':')	#-> ['chr1','1304646-1306645']
		chromosome = range[0]	#-> chr1
		range = range[1]	#-> 1304646-1306645
		range = range.split('-')
		range = map(int, range)	#-> [1304646, 1306645]
		prom_genome_start = range[0]	#-> 1304646
		prom_genome_end = range[1]	#-> 1306645
		strand = ls[5]	#-> 'strand=+'
		strand = strand.split('=')[1]	#-> '+'
		
		return acc, chromosome, prom_genome_start, prom_genome_end, strand
	
	def prom_seq_parse(self, curs, inputfile, output_table, organism, sequence_type, run_type=1):
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = fasta_block_iterator(inf)
		block_no = 0
		for block in iter:
			block = cStringIO.StringIO(block)
			header_line = block.readline()
			acc, chromosome, prom_genome_start, prom_genome_end, strand = self.parse_fasta_first_line(header_line)
			sequence = ''
			for line in block:
				sequence += line[:-1]
			row = [acc, chromosome, strand, prom_genome_start, prom_genome_end, organism, sequence_type, sequence]
			self.submit2prom_seq(curs, output_table, row)
			block_no += 1
			if self.report and block_no%500==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		
		del iter, inf
		sys.stderr.write("Done\n")
	
	def transfac_embl_parse(self, block, run_type):
		unit = transfac_embl_block()
		binding_region_list = []
		cell_source_list =[]
		comment_list = []
		expr_pattern_list = []
		external_database_links_list = []
		factors_list = []
		feature_list = []
		func_feature_list = []
		interacting_partners_list = []
		matrices_list = []
		method_list = []
		pwm_list = []
		reference_list = []
		sites_list = []
		struc_feature_list = []
		site_in_matrix_accs_list = []
		site_in_matrix_id_counter = 0
		block = cStringIO.StringIO(block)
		for line in block:
			try:
				line = line.replace("'", 'PRIME')	#replace the ' to avoid database-submit error
				line = line.replace('"'," DOUBLEPRIME ")
				if self.pwm_line_pattern.match(line):
					ls = line.split()
					pwm_list.append(map(float, ls[1:-1]))
					unit.consensus += ls[-1]
				if line[:4] == 'AC  ':
					unit.acc = line[4:-1]
				if line[:4] == 'ID  ':
					unit.id = line[4:-1]

				if line[:4] == 'BA  ':
					unit.basis += line[4:-1]
				if line[:4] == 'BC  ':
					unit.prom_classification = line[4:-1]
				if line[:4] == 'BF  ':
					factor_acc = line[4:-2].split('; ')[0]	#Watch: it's -2
					factors_list.append(factor_acc)
				if line[:4] == 'BR  ':
					fragment_acc = line[4:-2].split('; ')[0]	#Watch: it's -2
					binding_region_list.append(fragment_acc)
				if line[:4] == 'BS  ':
					if run_type ==3:	#matrix.dat has totally different format
						site_in_matrix_id_counter += 1
						site_in_matrix_acc = '%s_%s'%(unit.acc, site_in_matrix_id_counter)
						site_in_matrix_accs_list.append(site_in_matrix_acc)
						sites_seq, sites_acc, mt_start, len, gap_list,orientation = line[4:-2].replace('; ',';').split(';')	#replace '; ' with ';', then split with ';'
						if gap_list:
							gap_list = '{' + gap_list + '}'
						else:
							gap_list = '{}'
						unit.site_in_matrix_list.append([site_in_matrix_acc, unit.acc, sites_acc, sites_seq, \
							int(mt_start), int(len), gap_list, orientation])	#for table site_in_matrix
					else:
						sites_acc = line[4:-2].split('; ')[0]
					sites_list.append(sites_acc)
				if line[:4] == 'CC  ':
					comment_list.append(line[4:-1])
				if line[:4] == 'CD  ':
					unit.cell_description = line[4:-1]
				if line[:4] == 'CH  ':
					unit.chromosome = line[4:-1]
				if line[:4] == 'CL  ':
					if run_type==10:	#class.dat has different format
						unit.class_struc = line[4:-1]
					else:
						ls = line[4:-2].split('; ')
						if len(ls)>=3:
							unit.class_acc, unit.class_tree_id = ls[0],ls[2]	#no class_tree_id
						else:
							unit.class_acc = ls[0]
				if line[:4] == 'CN  ':
					unit.cell_negative = line[4:-1]
				if line[:4] == 'CP  ':
					unit.cell_positive = line[4:-1]
				if line[:4] == 'DE  ':
					unit.description = line[4:-1]
					if unit.description.find('Gene: ') != -1:
						unit.gene_acc = unit.description[:-1].split('Gene: ')[1]	#-1 is the period
				if line[:4] == 'DR  ':
					external_database_links_list.append(line[4:-2])
				if line[:4] == 'EL  ':
					unit.element = line[4:-1]
				if line[:4] == 'EX  ':
					expr_pattern_list.append(line[4:-1])
				if line[:4] == 'FA  ':
					if run_type ==2:	#in factor.dat, it's different
						unit.tf_name = line[4:-1]
					else:
						factor_acc = line[4:-2].split('; ')[0]	#Watch: it's -2
						factors_list.append(factor_acc)
				if line[:4] == 'FF  ':
					func_feature_list.append(line[4:-1])
				if line[:4] == 'FT  ':
					feature_list.append(line[4:-1])
				if line[:4] == 'GE  ':
					unit.gene_acc = line[4:-1].split('; ')[0]
				if line[:4] == 'HO  ':
					unit.homolog = line[4:-1]
				if line[:4] == 'IN  ':
					interacting_partners_acc = line[4:-1].split('; ')[0]
					interacting_partners_list.append(interacting_partners_acc)
				if line[:4] == 'MM  ':
					method_list.append(line[4:-1])
				if line[:4] == 'MX  ':
					matrix_acc = line[4:-1].split('; ')[0]
					matrices_list.append(matrix_acc)
				if line[:4] == 'NA  ':
					unit.tf_name = line[4:-1]
				if line[:4] == 'OS  ':
					organism = line[4:-1]
					organism = organism.split(', ')
					if len(organism)>1:	#only one name
						unit.organism = organism[1]
					else:
						unit.organism = organism[0]
				if line[:4] == 'RA  ':
					unit.ref_authors = line[4:-1]
				if line[:4] == 'RE  ':
					unit.gene_region = line[4:-1]
				if line[:4] == 'RG  ':
					unit.regulation = line[4:-1]
				if line[:4] == 'RL  ':
					unit.ref_journal = line[4:-1]
				if line[:4] == 'RN  ':
					ref_acc = line[4:-2].split()[1]	#Watch: it's line[4:-2]
					reference_list.append(ref_acc)

				if line[:4] == 'RT  ':
					unit.ref_title = line[4:-1]
				if line[:4] == 'RX  ':
					unit.ref_external_link = line[4:-2]

				if line[:4] == 'S1  ':
					unit.position_ref_point = line[4:-1]
				if line[:4] == 'SC  ':
					unit.sequence_source = line[4:-1]
				if line[:4] == 'SD  ':
					unit.short_description = line[4:-1]
				if line[:4] == 'SF  ':
					if run_type==2:
						struc_feature_list.append(line[4:-1])
					else:
						unit.position_start = int(line[4:-1])
				if line[:4] == 'SO  ':
					cell_source_acc = line[4:-1].split(';')[0]
					cell_source_list.append(cell_source_acc)
					unit.cell_factor_source = line[4:-1]
				if line[:4] == 'SQ  ':
					if line[-1]=='.':
						unit.sequence += line[4:-2]	#for site.dat, remove the trailing period
					else:
						unit.sequence += line[4:-1]

				if line[:4] == 'ST  ':
					unit.position_end = int(line[4:-1])
				if line[:4] == 'SY  ':
					unit.synonyms = line[4:-1]
				if line[:4] == 'SZ  ':
					unit.size = line[4:-1]
				if line[:4] == 'TY  ':
					unit.type = line[4:-1]
			
			except:
				print 'Except: %s'%repr(sys.exc_info()[0])
				print line
				sys.exit(2)		
		if binding_region_list:
			unit.binding_region = '{' + ','.join(binding_region_list) + '}'
		if cell_source_list:
			unit.cell_source = '{' + ','.join(cell_source_list) + '}'
		if comment_list:
			unit.comment = '||'.join(comment_list)
		if expr_pattern_list:
			unit.expr_pattern = '{"' + '","'.join(expr_pattern_list) + '"}'	#array '{}' could only contain ". only for varchar[]
		if external_database_links_list:
			unit.external_database_links = '{"' + '","'.join(external_database_links_list) + '"}'
		if factors_list:
			unit.factors = '{' + ','.join(factors_list) + '}'
		if feature_list:
			unit.feature = '{"' + '","'.join(feature_list) + '"}'
		if func_feature_list:
			unit.func_feature = '{"' + '","'.join(func_feature_list) + '"}'
		if interacting_partners_list:
			unit.interacting_partners = '{' + ','.join(interacting_partners_list) + '}'
		if matrices_list:
			unit.matrices = '{' + ','.join(matrices_list) + '}'
		if method_list:
			unit.method = '{"' + '","'.join(method_list) + '"}'
		if pwm_list:
			pwm = repr(pwm_list)
			pwm = pwm.replace('[','{')
			unit.pwm = pwm.replace(']','}')
		if reference_list:
			unit.reference = '{' + ','.join(reference_list) + '}'
		if site_in_matrix_accs_list:
			unit.site_in_matrix_accs = '{' + ','.join(site_in_matrix_accs_list) + '}'
		if sites_list:
			unit.sites = '{' + ','.join(sites_list) + '}'
		if struc_feature_list:
			unit.struc_feature = '{"' + '","'.join(struc_feature_list) + '"}'
		return unit
		
	def submit2factor(self, curs, table_name, row):
		curs.execute("insert into %s(tf_acc, tf_id, tf_name, tf_syn, organism, gene_acc, homolog,\
			class_acc, class_tree_id, size, sequence, sequence_source, feature, struc_feature,\
			cell_positive, cell_negative, expr_pattern, func_feature, interacting_partners, matrices,\
			sites, binding_region, external_database_links, reference) values('%s')"%(table_name, "', '".join(row)))
	
	def factor_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=2):
		"""
		09-15-05
			reference is a list of ref_acc
		"""
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = unigene_data_block_iterator(inf)
		block_no = 0
		for block in iter:
			block_no += 1
			if block_no==1:	#skip the first block
				continue
			
			unit = self.transfac_embl_parse(block, run_type)
			row = [unit.acc, unit.id, unit.tf_name, unit.synonyms, unit.organism, unit.gene_acc, unit.homolog,\
				unit.class_acc, unit.class_tree_id, unit.size, unit.sequence, unit.sequence_source, unit.feature,\
				unit.struc_feature, unit.cell_positive, unit.cell_negative, unit.expr_pattern, unit.func_feature,\
				unit.interacting_partners, unit.matrices, unit.sites, unit.binding_region, unit.external_database_links,\
				unit.reference]
			self.submit2factor(curs, output_table, row)
			
			if self.report and block_no%50==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		sys.stderr.write("Done\n")
	
	def submit2matrix(self, curs, table_name, row, site_in_matrix_list):
		"""
		09-16-05
			add site_in_matrix_accs to output_table('matrix')
			add site_in_matrix_list to submit to site_in_matrix
		"""
		curs.execute("insert into %s(mt_acc, mt_id, tf_name, tf_description, factors, pwm, consensus,\
			basis, sites, site_in_matrix_accs, comment, reference) values('%s')"%(table_name, "', '".join(row)))
		for row in site_in_matrix_list:
			curs.execute("insert into site_in_matrix(acc, mt_acc, site_acc, sequence, mt_start, len, \
			gap_list, orientation) values ('%s', %s, %s, '%s')"%("','".join(row[:4]), row[4], row[5], "','".join(row[6:])))
	
	def matrix_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=3):
		"""
		09-15-05
			reference is a list of ref_acc
			add pwm and consensus
		09-16-05
			add site_in_matrix_accs to output_table('matrix')
			add site_in_matrix_list to submit to site_in_matrix
		"""
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = unigene_data_block_iterator(inf)
		block_no = 0
		pwm_line_pattern = re.compile(r'\d\d  ')	#used to identify pwm lines, first two characters are numbers
		for block in iter:
			block_no += 1
			if block_no == 1:	#skip the first block
				continue
			unit = self.transfac_embl_parse(block, run_type)
			row = [unit.acc, unit.id, unit.tf_name, unit.description, unit.factors, unit.pwm, unit.consensus, \
				unit.basis, unit.sites, unit.site_in_matrix_accs, unit.comment, unit.reference]
			self.submit2matrix(curs, output_table, row, unit.site_in_matrix_list)
			
			if self.report and block_no%500==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		sys.stderr.write("Done\n")
	
	def submit2binding_site(self, curs, table_name, row):
		curs.execute("insert into %s(mt_id, strand, bs_disp_start, bs_disp_end, \
			prom_id, core_similarity_score, matrix_similarity_score, sequence)\
			values('%s', %s, '%s')"%(table_name, "', '".join(row[:2]), repr(row[2:-1])[1:-1], row[-1]))
	
	def binding_site_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=4):
		"""
		09-14-05
			some vacant values(chromosome, bs_genome_start, bs_genome_end...) will be
				filled in through database's trigger
		"""
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = match_block_iterator(inf)
		block_no = 0
		for block in iter:
			mt_id = ''
			strand = ''
			bs_disp_start = -1
			bs_disp_end = -1
			prom_id = None
			core_similarity_score = None
			matrix_similarity_score = None
			sequence = ''
			
			block = cStringIO.StringIO(block)
			line = block.readline()
			prom_id = int(line.split()[-1])	#\n is discarded automatically by split()
			for line in block:
				try:
					ls = line[:-1].split('|')
					mt_id = ls[0].replace(' ','')	#remove spaces
					bs_disp_start_strand = ls[1].replace(' ','')
					bs_disp_start = int(bs_disp_start_strand[:-3])
					strand = bs_disp_start_strand[-2]
					core_similarity_score = float(ls[2])
					matrix_similarity_score = float(ls[3])
					sequence = ls[4].replace(' ','')
					bs_disp_end = bs_disp_start + len(sequence) - 1
					row = [mt_id, strand, bs_disp_start, bs_disp_end, prom_id, \
						core_similarity_score, matrix_similarity_score, sequence]
					self.submit2binding_site(curs, output_table, row)
				except:
					print line
					print ls
					sys.exit(2)
			
			block_no += 1
			if self.report and block_no%1000==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		sys.stderr.write("Done\n")
	
	def submit2reference(self, curs, table_name, row):
		curs.execute("insert into %s(ref_acc, ref_external_link, ref_authors, ref_title, ref_journal)\
			values('%s')"%(table_name, "', '".join(row)))
	
	def reference_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=5):
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = unigene_data_block_iterator(inf)
		block_no = 0
		for block in iter:
			block_no += 1
			if block_no == 1:	#skip the first block
				continue
			ref_acc = ''
			ref_external_link = ''
			ref_authors = ''
			ref_title = ''
			ref_journal	= ''
			block = cStringIO.StringIO(block)
			for line in block:
				line = line.replace("'", 'PRIME')	#replace the ' to avoid database-submit error
				if line[:4] == 'AC  ':
					ref_acc = line[4:-1]
				if line[:4] == 'RX  ':
					ref_external_link = line[4:-1]
				if line[:4] == 'RA  ':
					ref_authors = line[4:-1]
				if line[:4] == 'RT  ':
					ref_title = line[4:-1]
				if line[:4] == 'RL  ':
					ref_journal = line[4:-1]
			row = [ref_acc, ref_external_link, ref_authors, ref_title, ref_journal]
			self.submit2reference(curs, output_table, row)
			
			if self.report and block_no%500==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		sys.stderr.write("Done\n")
	
	def submit2site(self, curs, table_name, row):
		curs.execute("insert into %s(acc, acc_id, type, description, gene_acc, organism, gene_region, \
			sequence, element, position_ref_point, position_start, position_end, cell_source, method, comment, \
			external_database_links, reference) values('%s', %s, %s, '%s')"%(table_name, "', '".join(row[:10]), \
			row[10], row[11], "', '".join(row[12:])))
	
	def site_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=6):
		"""
		09-15-05
		"""
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = unigene_data_block_iterator(inf)
		block_no = 0
		pwm_line_pattern = re.compile(r'\d\d  ')	#used to identify pwm lines, first two characters are numbers
		for block in iter:
			block_no += 1
			if block_no == 1:	#skip the first block
				continue
			
			unit = self.transfac_embl_parse(block, run_type)
			row = [unit.acc, unit.id, unit.type, unit.description, unit.gene_acc, unit.organism, unit.gene_region, unit.sequence,\
				unit.element, unit.position_ref_point, unit.position_start, unit.position_end, unit.cell_source, unit.method, \
				unit.comment, unit.external_database_links, unit.reference]
			self.submit2site(curs, output_table, row)
			
			if self.report and block_no%500==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		sys.stderr.write("Done\n")
	
	def submit2cell(self, curs, table_name, row):
		curs.execute("insert into %s(acc, cell_factor_source, organism, cell_description, sites, \
				binding_region, external_database_links) values('%s')"%(table_name, "', '".join(row)))
	
	def cell_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=7):
		"""
		09-15-05
		"""
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = unigene_data_block_iterator(inf)
		block_no = 0
		pwm_line_pattern = re.compile(r'\d\d  ')	#used to identify pwm lines, first two characters are numbers
		for block in iter:
			block_no += 1
			if block_no == 1:	#skip the first block
				continue
			
			unit = self.transfac_embl_parse(block, run_type)
			row = [unit.acc, unit.cell_factor_source, unit.organism, unit.cell_description, unit.sites, \
				unit.binding_region, unit.external_database_links]
			self.submit2cell(curs, output_table, row)
			
			if self.report and block_no%500==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		sys.stderr.write("Done\n")
	
	def submit2fragment(self, curs, table_name, row):
		curs.execute("insert into %s(acc, acc_id, description, gene_acc, organism, sequence, sequence_source, \
				factors, method, external_database_links, reference) values('%s')"%(table_name, "', '".join(row)))
	
	def fragment_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=8):
		"""
		09-15-05
		"""
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = unigene_data_block_iterator(inf)
		block_no = 0
		pwm_line_pattern = re.compile(r'\d\d  ')	#used to identify pwm lines, first two characters are numbers
		for block in iter:
			block_no += 1
			if block_no == 1:	#skip the first block
				continue
			
			unit = self.transfac_embl_parse(block, run_type)
			row = [unit.acc, unit.id, unit.description, unit.gene_acc, unit.organism, unit.sequence, unit.sequence_source, \
				unit.factors, unit.method, unit.external_database_links, unit.reference]
			self.submit2fragment(curs, output_table, row)
			
			if self.report and block_no%500==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		sys.stderr.write("Done\n")
	
	def submit2gene(self, curs, table_name, row):
		curs.execute("insert into %s(acc, acc_id, short_description, description, synonyms, organism,\
			chromosome, prom_classification, regulation, binding_region, factors, external_database_links, \
			reference) values('%s')"%(table_name, "', '".join(row)))
	
	def gene_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=9):
		"""
		09-15-05
		"""
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = unigene_data_block_iterator(inf)
		block_no = 0
		pwm_line_pattern = re.compile(r'\d\d  ')	#used to identify pwm lines, first two characters are numbers
		for block in iter:
			block_no += 1
			if block_no == 1:	#skip the first block
				continue
			
			unit = self.transfac_embl_parse(block, run_type)
			row = [unit.acc, unit.id, unit.short_description, unit.description, unit.synonyms,\
				unit.organism, unit.chromosome, unit.prom_classification, unit.regulation, unit.binding_region,\
				unit.factors, unit.external_database_links, unit.reference]
			self.submit2gene(curs, output_table, row)
			
			if self.report and block_no%500==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		sys.stderr.write("Done\n")
	
	def submit2class(self, curs, table_name, row):
		curs.execute("insert into %s(acc, acc_id, class_struc, description, factors, external_database_links, \
			reference) values('%s')"%(table_name, "', '".join(row)))
	
	def class_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=10):
		"""
		09-15-05
		"""
		sys.stderr.write("Parsing for %s...\n"%output_table)
		inf = open(inputfile,'r')
		iter = unigene_data_block_iterator(inf)
		block_no = 0
		pwm_line_pattern = re.compile(r'\d\d  ')	#used to identify pwm lines, first two characters are numbers
		for block in iter:
			block_no += 1
			if block_no == 1:	#skip the first block
				continue
			
			unit = self.transfac_embl_parse(block, run_type)
			row = [unit.acc, unit.id, unit.class_struc, unit.comment, \
				unit.factors, unit.external_database_links, unit.reference]	#comment is description
			self.submit2class(curs, output_table, row)
			
			if self.report and block_no%500==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		sys.stderr.write("Done\n")
		
	def binding_site_easy_parse(self, curs, inputfile, output_table, organism_given, sequence_type, run_type=11):
		"""
		01-03-06
		an easy version to deal with output of TopMatchResult.py
		"""
		sys.stderr.write("Easy binding_site parsing for %s...\n"%output_table)
		reader = csv.reader(open(inputfile,'r'), delimiter='\t')
		counter = 0
		for row in reader:
			try:
				mt_id, prom_id, strand, bs_disp_start, bs_disp_end, core_similarity_score, \
					matrix_similarity_score, sequence = row
				new_row = [mt_id, strand, bs_disp_start, bs_disp_end, prom_id, \
					core_similarity_score, matrix_similarity_score, sequence]
				self.submit2binding_site(curs, output_table, new_row)
			except:
				print line
				print ls
				sys.exit(2)
			
			counter += 1
			if self.report and counter%5000==0:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		del reader
		sys.stderr.write("Done\n")
	
	def run(self):
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		self.parser_dict[self.type](curs, self.inputfile, self.output_table_dict[self.type], self.organism, self.sequence_type)
		if self.commit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:p:g:y:brch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'transfac'
	inputfile = None
	output_table = None
	sequence_type = 1
	organism = 'hs'
	type = 1
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
		elif opt in ("-i",):
			inputfile = arg
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-p",):
			sequence_type = int(arg)
		elif opt in ("-g",):
			organism = arg
		elif opt in ("-y",):
			type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1

	if inputfile and hostname and dbname and schema:
		instance = transfacdb(hostname, dbname, schema, inputfile, output_table, \
			sequence_type, organism, type, debug, report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
