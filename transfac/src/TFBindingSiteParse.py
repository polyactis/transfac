#!/usr/bin/env python
"""
Usage: TFBindingSiteParse.py [OPTIONS] -k SCHEMA -i INPUT_FILE -y PARSER_TYPE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...	INPUT_FILE
	-g ...	the two-letter abbreviation of organism, hs(default)
		only for type =1
	-y ...	parser type(1, default)
	-b, --debug	just test running the program, no daytime restriction
	-r, --report	report the progress (time before each query)
	-c, --commit	commit this database transaction
	-h, --help	show this help

Examples:
	TFBindingSiteParse.py -k harbison2004 -i 
		/usr/local/research_data/Harbison2004/p_value_excel/pvalbygene_forpaper_abbr_YPD.csv
		-g sc -c  >/tmp/yuhuang/harbison2004_missed_orf
		
	~/script/transfac/src/TFBindingSiteParse.py -k sgd_regulatory -i ./scerevisiae_regulatory.gff -g sc -c -y 3

Description:
	Program to construct tfbs.sql from inputfiles(from paper, other database)

	parser type
		1: harbison2004, 2:cisred.features, 3: sgd_regulatory, 4: ucsc_tfbs_conserved_parse
		info for 4: http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=68793971&c=chrX&g=tfbsConsSites
"""

import sys, getopt, os, csv, cStringIO, re
sys.path += [os.path.join(os.path.expanduser('~/script'))]
from annot.bin.codense.common import db_connect, org_short2long, org2tax_id, \
	get_gene_symbol2gene_id, get_ensembl_id2gene_id, get_entrezgene_coord,\
	get_entrezgene_annotated_anchor, get_sequence_segment, calculate_rome_number
from sets import Set
from Bio.Seq import Seq	#12-18-05 for reverse_complement()
from pymodule import PassingData

class binding_site_attribute:
	def __init__(self, mt_id='', chromosome='', strand='', bs_genome_start=-1, bs_genome_end=-1,\
		bs_disp_start=-1, bs_disp_end=-1, prom_id=-1, tax_id=-1, core_similarity_score=-1, matrix_similarity_score=-1,
		sequence='', comment=''):
		self.mt_id = mt_id
		self.chromosome = chromosome
		self.strand = strand
		self.bs_genome_start = bs_genome_start
		self.bs_genome_end = bs_genome_end
		self.bs_disp_start = bs_disp_start
		self.bs_disp_end = bs_disp_end
		self.prom_id = prom_id
		self.tax_id = tax_id
		self.core_similarity_score = core_similarity_score
		self.matrix_similarity_score  = matrix_similarity_score
		self.sequence = sequence
		self.comment = comment

class TFBindingSiteParse(object):
	"""
	12-18-05
		add sgd_regulatory related functions
	"""
	def __init__(self, hostname='zhoudb', dbname='mdb', schema='', inputfile=None, \
		organism='hs', type=1, debug=0, report=0, commit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputfile = inputfile
		self.organism = organism
		self.type = int(type)
		self.debug = int(debug)
		self.report = int(report)
		self.commit = int(commit)
		
		self.tax_id = org2tax_id(org_short2long(self.organism))
		self.parser_dict = {1: self.harbison2004_parse,
			2: self.cisred_parse,
			3: self.sgd_regulatory_parse,
			4: self.ucsc_tfbs_conserved_parse}
	
	def submit2binding_site(self, curs, binding_site_attr, table_name='binding_site'):
		"""
		12-10-05
			add bs_genome_start and bs_genome_end
		"""
		curs.execute("insert into %s(mt_id, chromosome, strand, bs_genome_start, bs_genome_end, \
			bs_disp_start, bs_disp_end, \
			prom_id, tax_id, core_similarity_score, matrix_similarity_score, sequence, comment)\
			values('%s', '%s', '%s', %s, %s, \
			%s, %s,\
			%s, %s, %s, %s, '%s', '%s')"%(table_name, binding_site_attr.mt_id, binding_site_attr.chromosome, \
			binding_site_attr.strand,\
			binding_site_attr.bs_genome_start, binding_site_attr.bs_genome_end,\
			binding_site_attr.bs_disp_start, binding_site_attr.bs_disp_end, binding_site_attr.prom_id,\
			binding_site_attr.tax_id, binding_site_attr.core_similarity_score, \
			binding_site_attr.matrix_similarity_score, binding_site_attr.sequence, binding_site_attr.comment))
	
	def expand_matrix_table(self, curs, mt_id_set, tax_id, matrix_table='matrix'):
		sys.stderr.write("\t Expanding matrix table")
		curs.execute("select mt_id from matrix where tax_id=tax_id")
		rows = curs.fetchall()
		for row in rows:
			mt_id = row[0]
			if mt_id in mt_id_set:	#already present, remove it
				mt_id_set.remove(mt_id)
		#insert the remaining new mt_id into matrix
		for mt_id in mt_id_set:
			curs.execute("insert into %s(mt_id, tax_id) values('%s', %s)"%\
				(matrix_table, mt_id, tax_id))
		sys.stderr.write(".\n")
	
	def harbison2004_parse_header(self, row):
		sys.stderr.write("\t Parsing harbison2004 header")
		index2mt_id_condition = {}
		mt_id_set = Set()
		for i in range(3, len(row)):	#the 1st 3 fields are for target genes
			mt_id, condition = row[i].split('_')
			index2mt_id_condition[i] = [mt_id, condition]
			mt_id_set.add(mt_id)
		sys.stderr.write(".\n")
		return index2mt_id_condition, mt_id_set
	
	def harbison2004_parse(self, curs, inputfile, gene_symbol2gene_id, tax_id):
		"""
		12-09-05
			
			--harbison2004_parse_header()
			--expand_matrix_table()
			--submit2binding_site()
		"""
		sys.stderr.write("Parsing harbison2004 ...\n")
		inf = csv.reader(open(inputfile, 'r'), delimiter = '\t')
		#skip the 1st two lines
		inf.next()
		index2mt_id_condition, mt_id_set = self.harbison2004_parse_header(inf.next())
		self.expand_matrix_table(curs, mt_id_set, tax_id)
		
		for row in inf:
			gene_symbol, alias, gene_desc = row[:3]
			for i in range(3, len(row)):
				try:
					if row[i] and row[i]!='NaN':
						p_value = float(row[i])
						if p_value>0.001:	#too high p_value
							continue
						if gene_symbol not in gene_symbol2gene_id:
							sys.stdout.write("\t Warning: %s not in gene_symbol2gene_id.\n"%gene_symbol)
							continue
						binding_site_attr = binding_site_attribute(tax_id=tax_id)
						binding_site_attr.prom_id = gene_symbol2gene_id[gene_symbol]
						binding_site_attr.core_similarity_score = p_value
						binding_site_attr.mt_id = index2mt_id_condition[i][0]
						binding_site_attr.comment = index2mt_id_condition[i][1]
						self.submit2binding_site(curs, binding_site_attr)
				except:
					sys.stderr.write("Error while parsing: %s.\n"%repr(sys.exc_info()))
					sys.stderr.write("gene_symbol: %s, alias: %s, i: %s, row[i]: %s\n"%(gene_symbol, alias, i, row[i]))
					sys.exit(3)
		sys.stderr.write("Done.\n")
	
	def find_binding_site_disp_coord(self, binding_site_attr, entrezgene_strand, entrezgene_start, entrezgene_stop):
		"""
		12-10-05
		"""
		if entrezgene_strand == '+':
			binding_site_attr.bs_disp_start = binding_site_attr.bs_genome_start - entrezgene_start
			binding_site_attr.bs_disp_end = binding_site_attr.bs_genome_end - entrezgene_start
		else:
			binding_site_attr.bs_disp_start = entrezgene_stop - binding_site_attr.bs_genome_end
			binding_site_attr.bs_disp_end = entrezgene_stop - binding_site_attr.bs_genome_start
	
	def submit2matrix(self, curs, mt_id, tax_id, table_name='matrix'):
		curs.execute("insert into %s (mt_id, tax_id) values ('%s', %s)"%(table_name, mt_id, tax_id))
	
	def cisred_parse(self, curs, inputfile, ensembl_id2gene_id, tax_id):
		"""
		12-10-05
			
			--get_entrezgene_coord()
			--find_binding_site_disp_coord()
			--submit2matrix()
			--submit2binding_site()
		"""
		sys.stderr.write("Parsing cisred ...\n")
		inf = csv.reader(open(inputfile, 'r'), delimiter = '\t')
		#skip the 1st line
		inf.next()
		for row in inf:
			id, batch_id, seqname, source, feature, start, end, score, mi_score,\
				strand, frame, ensembl_gene_id, consensus, reverse_consensus = row
			if ensembl_gene_id not in ensembl_id2gene_id:
				sys.stdout.write("\t Warning: %s not in ensembl_id2gene_id.\n"%ensembl_gene_id)
				continue
			mt_id_not_submited = 1
			for gene_id in ensembl_id2gene_id[ensembl_gene_id]:
				gene_id_coord = get_entrezgene_coord(curs, gene_id)
				if gene_id_coord is None:
					sys.stdout.write("\t\t Warning: %s not in entrezgene_mapping.\n"%gene_id)
					continue
				genomic_gi, chromosome, entrezgene_strand, entrezgene_start, entrezgene_stop = gene_id_coord
				if chromosome != seqname:
					sys.stdout.write("\t\t Warning: entrezgene %s, chromosome %s diffs from ensembl_gene_id: %s chromosome %s.\n"%\
						(gene_id, chromosome, ensembl_gene_id, seqname))
					continue
				binding_site_attr = binding_site_attribute(tax_id=tax_id)
				binding_site_attr.chromosome = chromosome
				binding_site_attr.strand = strand
				binding_site_attr.mt_id = 'craHsap%s'%id
				if mt_id_not_submited:
					self.submit2matrix(curs, binding_site_attr.mt_id, tax_id)	#two gene_id for one ensembl_gene_id
					mt_id_not_submited = 0
				binding_site_attr.bs_genome_start = int(start)
				binding_site_attr.bs_genome_end = int(end)
				binding_site_attr.prom_id = gene_id
				binding_site_attr.core_similarity_score = float(score)
				binding_site_attr.matrix_similarity_score = float(mi_score)
				binding_site_attr.comment = feature
				if consensus !='NULL':
					binding_site_attr.sequence = consensus
				elif reverse_consensus != 'NULL':
					binding_site_attr.sequence = reverse_consensus
				self.find_binding_site_disp_coord(binding_site_attr, entrezgene_strand, entrezgene_start, entrezgene_stop)
				self.submit2binding_site(curs, binding_site_attr)
		sys.stderr.write("Done.\n")
	
	def parse_sgd_regulatory_attribute(self, attribute):
		"""
		12-17-05
			to get mt_id and dbxref
			i.e.
	ID=SWI6-binding-site-S000086008;Name=SWI6-binding-site-S000086008;Note=SWI6%20binding%20site;dbxref=SGD:S000086008
		"""
		row = attribute.split(';')
		for entry in row:
			name, desc = entry.split('=')
			if name=='ID':
				end_index = desc.find('-binding-site')
				mt_id = desc[:end_index]
		dbxref = row[-1]	#the last one
		return mt_id, dbxref
	
	def addTargetGeneToGeneLs(cls, bs_anchor, gene_anchor, gene_id, target_gene_ls, left_or_right, \
							upstream_or_downstream=None, max_upstream_distance=-1, max_downstream_distance=-1):
		"""
		2008-08-19
			max_upstream_distance>-1 => max_upstream_distance>=0. more precise conditions.
		2008-08-13
			if max_upstream_distance or max_downstream_distance >0, it will take all genes within that distance, rather than the closest
		"""
		good_target_gene = 0
		abs_distance = abs(bs_anchor-gene_anchor)
		if upstream_or_downstream =='upstream':
			disp_pos = -abs_distance
		elif upstream_or_downstream=='downstream':
			disp_pos = +abs_distance
		else:
			disp_pos = abs_distance
		if upstream_or_downstream=='upstream' and max_upstream_distance>=0:
			if abs_distance<=max_upstream_distance:
				good_target_gene = 1
		elif upstream_or_downstream=='downstream' and max_downstream_distance>=0:
			if abs_distance<=max_downstream_distance:
				good_target_gene = 1
		elif target_gene_ls:	#code below is to take whatever closest to the bs_anchor
			if gene_anchor==target_gene_ls[0][0]:	#as close as the other
				good_target_gene = 1
			elif left_or_right=='right' and gene_anchor>target_gene_ls[0][0]:	#closer than the incumbent, replace it
				good_target_gene = 2
			elif left_or_right=='left' and gene_anchor<target_gene_ls[0][0]:	#closer than the incumbent, replace it
				good_target_gene = 2
		else:	#nothing in target_gene_ls and the two maximum distance isn't used, take whatever
			good_target_gene = 1
		
		if not upstream_or_downstream:	#copy left_or_right if it's nothing
			upstream_or_downstream = left_or_right
		if good_target_gene == 1:
				target_gene_ls.append((gene_anchor, gene_id, left_or_right, upstream_or_downstream, disp_pos))
		elif good_target_gene==2:
			target_gene_ls = [(gene_anchor, gene_id, left_or_right, upstream_or_downstream, disp_pos)]
	
	addTargetGeneToGeneLs = classmethod(addTargetGeneToGeneLs)
	
	def return_target_gene_ls(cls, regulatory_coord, chromosome2anchor_gene_tuple_ls, gene_id2coord, max_upstream_distance=-1, max_downstream_distance=-1):
		"""
		2008-08-13
			if max_upstream_distance or max_downstream_distance >0, it will take all genes within that distance, rather than the closest
		2008-08-12
			chromosome2anchor_gene_tuple_ls only contains a unique gene once
			reorganize the way to find target genes from which regulatory_coord is either upstream or downstream or left or right.
			in each category, only closest gene is added.
		12-17-05
			very compicated issue
			draw plots to understand it
		"""
		chromosome, bs_genome_start, bs_genome_end = regulatory_coord
		anchor_gene_tuple_ls = chromosome2anchor_gene_tuple_ls[chromosome]
		regulatory_is_right_upstream_target_gene_ls = []	#
		regulatory_is_right_downstream_target_gene_ls = []
		regulatory_is_right_target_gene_ls = []
		regulatory_is_left_upstream_target_gene_ls = []
		regulatory_is_left_downstream_target_gene_ls = []
		regulatory_is_left_target_gene_ls = []
		regulatory_touch_target_gene_ls = []	#might have several genes covering the regulatory sequence
		
		for anchor_gene_tuple in anchor_gene_tuple_ls:
			anchor, gene_id = anchor_gene_tuple	#anchor is useless
			gene_start, gene_stop, gene_strand, gene_genomic_gi = gene_id2coord[gene_id]
			if gene_strand == '-':
				if gene_stop<bs_genome_start:	#before
					cls.addTargetGeneToGeneLs(bs_genome_start, gene_stop, gene_id, regulatory_is_right_upstream_target_gene_ls, \
											'right', 'upstream', max_upstream_distance, max_downstream_distance)
				elif gene_stop>=bs_genome_start and gene_start<=bs_genome_end:
					regulatory_touch_target_gene_ls.append((gene_stop, gene_id, 'touch', 'touch', gene_stop-bs_genome_start))
				elif gene_start>bs_genome_end:
					cls.addTargetGeneToGeneLs(bs_genome_end, gene_start, gene_id, regulatory_is_left_downstream_target_gene_ls, \
											'left', 'downstream', max_upstream_distance, max_downstream_distance)
			elif gene_strand=='+':
				if gene_start>bs_genome_end:
					cls.addTargetGeneToGeneLs(bs_genome_end, gene_start, gene_id, regulatory_is_left_upstream_target_gene_ls, \
											'left', 'upstream', max_upstream_distance, max_downstream_distance)
				elif gene_start<=bs_genome_end and gene_stop>=bs_genome_start:
					regulatory_touch_target_gene_ls.append((gene_start, gene_id, 'touch', 'touch', bs_genome_start-gene_start))
				elif gene_stop<bs_genome_start:
					cls.addTargetGeneToGeneLs(bs_genome_start, gene_stop, gene_id, regulatory_is_right_downstream_target_gene_ls, \
											'right', 'downstream', max_upstream_distance, max_downstream_distance)
			else:	#no strand info
				if gene_stop<bs_genome_start:
					cls.addTargetGeneToGeneLs(bs_genome_start, gene_stop, gene_id, regulatory_is_right_target_gene_ls, \
											'right', '', max_upstream_distance, max_downstream_distance)
				elif gene_start<=bs_genome_end and gene_stop>=bs_genome_start:
					regulatory_touch_target_gene_ls.append((gene_start, gene_id, 'touch', 'touch', bs_genome_start-gene_start))
				elif gene_start>bs_genome_end:
					cls.addTargetGeneToGeneLs(bs_genome_end, gene_start, gene_id, regulatory_is_left_target_gene_ls, \
											'left', '', max_upstream_distance, max_downstream_distance)
		
		pdata = PassingData(regulatory_is_right_upstream_target_gene_ls = regulatory_is_right_upstream_target_gene_ls,	\
						regulatory_is_right_downstream_target_gene_ls = regulatory_is_right_downstream_target_gene_ls,\
						regulatory_is_right_target_gene_ls = regulatory_is_right_target_gene_ls,
						regulatory_is_left_upstream_target_gene_ls = regulatory_is_left_upstream_target_gene_ls,
						regulatory_is_left_downstream_target_gene_ls = regulatory_is_left_downstream_target_gene_ls, 
						regulatory_is_left_target_gene_ls = regulatory_is_left_target_gene_ls,
						regulatory_touch_target_gene_ls = regulatory_touch_target_gene_ls)
		
		return pdata
	
	return_target_gene_ls = classmethod(return_target_gene_ls)
	
	def sgd_regulatory_parse(self, curs, inputfile, gene_symbol2gene_id, tax_id):
		"""
		2008-08-12
			return_target_gene_ls() returns a different data structure.
		12-17-05
			
			--get_entrezgene_annotated_anchor()
			--calculate_rome_number()
			
			--parse_sgd_regulatory_attribute()
			--submit2matrix()
			--return_target_gene_ls()
			--find_binding_site_disp_coord()
			--get_sequence_segment()
			--submit2binding_site()
		"""
		sys.stderr.write("Parsing sgd_regulatory ...\n")
		chromosome2anchor_gene_tuple_ls, gene_id2coord = get_entrezgene_annotated_anchor(curs, tax_id)
		inf = open(inputfile, 'r')
		mt_id_set = Set()
		line_no = 0
		for line in inf:
			line_no+=1
			if line[0]=='#':	#skip the comment lines
				continue
			row = line[:-1].split('\t')
			try:
				seqname, source, feature, start, end, score, strand, frame, attribute = row
				chromosome_rome_number = seqname[3:]
				if feature=='TF_binding_site':
					binding_site_attr = binding_site_attribute(tax_id=tax_id)
					chromosome = calculate_rome_number(chromosome_rome_number)
					if chromosome:
						binding_site_attr.chromosome = repr(chromosome)	#varchar type in database
					else:	#it might be something like Mito, return None
						binding_site_attr.chromosome = chromosome_rome_number
					binding_site_attr.strand = strand
					mt_id, dbxref = self.parse_sgd_regulatory_attribute(attribute)
					if mt_id not in mt_id_set:
						self.submit2matrix(curs, mt_id, tax_id)
						mt_id_set.add(mt_id)
					binding_site_attr.mt_id = mt_id
					binding_site_attr.bs_genome_start = int(start)
					binding_site_attr.bs_genome_end = int(end)
					
					regulatory_coord = (binding_site_attr.chromosome, \
						binding_site_attr.bs_genome_start, binding_site_attr.bs_genome_end)
					pdata = self.return_target_gene_ls(regulatory_coord, chromosome2anchor_gene_tuple_ls, gene_id2coord)
					if pdata.regulatory_touch_target_gene_ls:
						target_gene_ls_type = 'touch'
						target_gene_ls = pdata.regulatory_touch_target_gene_ls
					elif pdata.regulatory_is_left_upstream_target_gene_ls or pdata.regulatory_is_right_upstream_target_gene_ls:
						target_gene_ls_type = 'upstream'
						target_gene_ls = pdata.regulatory_is_left_upstream_target_gene_ls + pdata.regulatory_is_right_upstream_target_gene_ls
					
					binding_site_attr.comment = '%s;%s'%(dbxref, target_gene_ls_type)
					for target_gene_tuple in target_gene_ls:
						anchor, gene_id = target_gene_tuple[:2]
						binding_site_attr.prom_id = gene_id
						gene_start, gene_stop, gene_strand, gene_genomic_gi = gene_id2coord[gene_id]
						self.find_binding_site_disp_coord(binding_site_attr, gene_strand, gene_start, gene_stop)
						binding_site_attr.sequence = get_sequence_segment(curs, gene_genomic_gi, \
							binding_site_attr.bs_genome_start, binding_site_attr.bs_genome_end)
						if strand == '-':	#reverse_complement
							seq = Seq(binding_site_attr.sequence)
							binding_site_attr.sequence = seq.reverse_complement().tostring()
						self.submit2binding_site(curs, binding_site_attr)
			except:
				print "line_no:", line_no
				print line
		sys.stderr.write("Done.\n")
	
	def ucsc_tfbs_conserved_parse(self, curs, inputfile, gene_symbol2gene_id, tax_id):
		"""
		2008-08-12
			return_target_gene_ls() returns a different data structure.
		03-29-06
			
			--get_entrezgene_annotated_anchor()
			
			--submit2matrix()
			--return_target_gene_ls()
			--find_binding_site_disp_coord()
			--get_sequence_segment()
			--submit2binding_site()
		"""
		sys.stderr.write("Parsing ucsc_tfbs_conserved ...\n")
		chromosome2anchor_gene_tuple_ls, gene_id2coord = get_entrezgene_annotated_anchor(curs, tax_id)
		inf = open(inputfile, 'r')
		mt_id_set = Set()
		line_no = 0
		for line in inf:
			line_no+=1
			if line[0]=='#':	#skip the comment lines
				continue
			row = line[:-1].split('\t')
			try:
				bin, chromosome, start, end, name, score, strand, zscore = row
				chromosome = chromosome[3:]
				binding_site_attr = binding_site_attribute(tax_id=tax_id)
				if chromosome.find('random')==-1:	#03-29-06 not random sequences
					binding_site_attr.chromosome = chromosome
					binding_site_attr.strand = strand
					mt_id = name
					if mt_id not in mt_id_set:
						self.submit2matrix(curs, mt_id, tax_id)
						mt_id_set.add(mt_id)
					binding_site_attr.mt_id = mt_id
					binding_site_attr.bs_genome_start = int(start)
					binding_site_attr.bs_genome_end = int(end)
					binding_site_attr.matrix_similarity_score = float(score)
					binding_site_attr.core_similarity_score = float(zscore)
					
					regulatory_coord = (binding_site_attr.chromosome, \
						binding_site_attr.bs_genome_start, binding_site_attr.bs_genome_end)
					target_gene_ls, target_gene_ls_type = self.return_target_gene_ls(regulatory_coord, \
						chromosome2anchor_gene_tuple_ls, gene_id2coord)
					pdata = self.return_target_gene_ls(regulatory_coord, chromosome2anchor_gene_tuple_ls, gene_id2coord)
					if pdata.regulatory_touch_target_gene_ls:
						target_gene_ls_type = 'touch'
						target_gene_ls = pdata.regulatory_touch_target_gene_ls
					elif pdata.regulatory_is_left_upstream_target_gene_ls or pdata.regulatory_is_right_upstream_target_gene_ls:
						target_gene_ls_type = 'upstream'
						target_gene_ls = pdata.regulatory_is_left_upstream_target_gene_ls + pdata.regulatory_is_right_upstream_target_gene_ls
					
					binding_site_attr.comment = 'bin:%s. %s'%(bin, target_gene_ls_type)
					for target_gene_tuple in target_gene_ls:
						anchor, gene_id = target_gene_tuple[:2]
						binding_site_attr.prom_id = gene_id
						gene_start, gene_stop, gene_strand, gene_genomic_gi = gene_id2coord[gene_id]
						self.find_binding_site_disp_coord(binding_site_attr, gene_strand, gene_start, gene_stop)
						binding_site_attr.sequence = get_sequence_segment(curs, gene_genomic_gi, \
							binding_site_attr.bs_genome_start, binding_site_attr.bs_genome_end)
						if strand == '-':	#reverse_complement
							seq = Seq(binding_site_attr.sequence)
							binding_site_attr.sequence = seq.reverse_complement().tostring()
						self.submit2binding_site(curs, binding_site_attr)
			except:
				print "line_no:", line_no
				print line
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		12-11-05
			
			--db_connect()
			--get_gene_symbol2gene_id()
			--get_ensembl_id2gene_id()
			--parser_dict[]
				--harbison2004_parse()
				--cisred_parse()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		if self.type==1:
			mapping_dict = get_gene_symbol2gene_id(curs, self.tax_id)
		elif self.type==2:
			mapping_dict = get_ensembl_id2gene_id(curs, self.tax_id)
		else:
			mapping_dict = None
		self.parser_dict[self.type](curs, self.inputfile, mapping_dict, self.tax_id)
		if self.commit:
			curs.execute("end")



if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:g:y:brch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = None
	inputfile = None
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
		elif opt in ("-i"):
			inputfile = arg
		elif opt in ("-g"):
			organism = arg
		elif opt in ("-y"):
			type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1

	if inputfile and hostname and dbname and schema:
		instance = TFBindingSiteParse(hostname, dbname, schema, inputfile,\
			organism, type, debug, report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
