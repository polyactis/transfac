#!/usr/bin/env python
"""

Examples:
	#postgresql
	GeneASNXML2gene_mapping.py -c -i *.xgs
	
	#package cElementTree is only installed in python 2.4
	python2.4 ./GeneASNXML2gene_mapping.py -z localhost -k genome -c -i /usr/local/research_data/ncbi/gene_2008_07_06/at.xgs -r

	# 2010-12-15
	GeneASNXML2gene_mapping.py -z localhost -k genome -d vervetdb  -i Pan_troglodytes.xgs -r -u yh -c

	#mysql 
	python2.4 GeneASNXML2gene_mapping.py -v mysql -u yh -z papaya -d genome -k "" -c -i at.xgs -r

Description:
	Program to parse the output of gene2xml and get location information
	for all genes and dump them into database.
	
	More info see "gene2xml info" and "Entrezgene XML tree hierarchy"
		of http://hto-pc44.usc.edu/research/annot/genomedb
	
"""

import sys, getopt, os
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
import xml.etree.cElementTree as ElementTree
from GenomeDB import GenomeDatabase, Gene, SequenceType, EntrezgeneType, \
	GeneSegment, GeneCommentaryType, GeneCommentary, AnnotAssembly
from datetime import datetime
from pymodule import PassingData
from pymodule.utils import returnAnyValueIfNothing
import traceback

class entrezgene_mapping_attr:
	"""
	2008-07-27
		deprecated
	"""
	def __init__(self):
		self.gene_id = None
		self.tax_id = None
		self.genomic_acc_ver = None
		self.genomic_gi = None
		self.strand = 'null'	#strand could be missing, +1 or -1 means plus or minus
		self.start = None
		self.stop = None
		self.mrna_acc_ver = None
		self.mrna_gi = None
		self.mrna_start = []
		self.mrna_stop = []
		self.cds_acc_ver = None
		self.cds_gi = None
		self.cds_start = []
		self.cds_stop = []
		self.comment = None

class GeneASNXML2gene_mapping:
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['graphdb', 'd', 1, 'database name', ],\
							('schema', 1, ): ['genome', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('inputfiles', 1, ):[None, 'i', 1, 'comma-separated input filenames, or a list of files'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-07-29
			use option_default_dict
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		if type(self.inputfiles)==str:
			self.inputfiles = self.inputfiles.split(',')
		
		self.find_info_dict = {'peptide':self.find_peptide_info,
			'rRNA':self.find_rna_info,
			'tRNA':self.find_rna_info,
			'snoRNA':self.find_rna_info,
			'miscRNA':self.find_rna_info,
			'snRNA':self.find_rna_info,
			'mRNA':self.find_mrna_info}
		
	def submit_to_entrezgene_mapping_table(self, curs, table, entrezgene_mapping):
		"""
		2008-07-27
			deprecated
		"""
		if entrezgene_mapping.mrna_acc_ver and entrezgene_mapping.cds_acc_ver:	#both rna and peptide present
			curs.execute("insert into %s(gene_id, tax_id, genomic_acc_ver, genomic_gi, strand, \
			start, stop, mrna_acc_ver, mrna_gi, mrna_start,\
			mrna_stop, cds_acc_ver, cds_gi, cds_start, cds_stop) values('%s', %s, '%s', %s, %s,\
			%s, %s, '%s', %s, ARRAY%s, \
			ARRAY%s, '%s', %s, ARRAY%s, ARRAY%s)"%(table,\
			entrezgene_mapping.gene_id, entrezgene_mapping.tax_id, entrezgene_mapping.genomic_acc_ver, entrezgene_mapping.genomic_gi, entrezgene_mapping.strand,\
			entrezgene_mapping.start, entrezgene_mapping.stop, entrezgene_mapping.mrna_acc_ver, entrezgene_mapping.mrna_gi, repr(entrezgene_mapping.mrna_start),\
			repr(entrezgene_mapping.mrna_stop), entrezgene_mapping.cds_acc_ver, entrezgene_mapping.cds_gi, repr(entrezgene_mapping.cds_start), repr(entrezgene_mapping.cds_stop)))
		elif entrezgene_mapping.mrna_acc_ver:	#only rna
			curs.execute("insert into %s(gene_id, tax_id, genomic_acc_ver, genomic_gi, strand, \
			start, stop, mrna_acc_ver, mrna_gi, mrna_start,\
			mrna_stop) values('%s', %s, '%s', %s, %s,\
			%s, %s, '%s', %s, ARRAY%s, \
			ARRAY%s)"%(table,\
			entrezgene_mapping.gene_id, entrezgene_mapping.tax_id, entrezgene_mapping.genomic_acc_ver, entrezgene_mapping.genomic_gi, entrezgene_mapping.strand,\
			entrezgene_mapping.start, entrezgene_mapping.stop, entrezgene_mapping.mrna_acc_ver, entrezgene_mapping.mrna_gi, repr(entrezgene_mapping.mrna_start),\
			repr(entrezgene_mapping.mrna_stop) ))
		elif entrezgene_mapping.cds_acc_ver:	#only peptide
			curs.execute("insert into %s(gene_id, tax_id, genomic_acc_ver, genomic_gi, strand, \
			start, stop, cds_acc_ver, cds_gi, cds_start, cds_stop) values('%s', %s, '%s', %s, %s,\
			%s, %s, '%s', %s, ARRAY%s, ARRAY%s)"%(table,\
			entrezgene_mapping.gene_id, entrezgene_mapping.tax_id, entrezgene_mapping.genomic_acc_ver, entrezgene_mapping.genomic_gi, entrezgene_mapping.strand,\
			entrezgene_mapping.start, entrezgene_mapping.stop, entrezgene_mapping.cds_acc_ver, entrezgene_mapping.cds_gi, repr(entrezgene_mapping.cds_start), repr(entrezgene_mapping.cds_stop)))
		else:	#only genomic location
			curs.execute("insert into %s(gene_id, tax_id, genomic_acc_ver, genomic_gi, strand, \
			start, stop) values('%s', %s, '%s', %s, %s,\
			%s, %s)"%(table,\
			entrezgene_mapping.gene_id, entrezgene_mapping.tax_id, entrezgene_mapping.genomic_acc_ver, entrezgene_mapping.genomic_gi, entrezgene_mapping.strand,\
			entrezgene_mapping.start, entrezgene_mapping.stop))
	
	def is_gi_valid_in_annot_assembly_table(self, curs, gi, annot_assembly_table):
		"""
		2008-07-28
			deprecated
		"""
		curs.execute("select gi from %s where gi=%s"%(annot_assembly_table, gi))
		rows = curs.fetchall()
		if rows:
			return 1
		else:
			return 0
	
	def returnSeqLocInt(self, elem):
		"""
		2008-07-29
		"""
		start_index = elem.findtext('Seq-interval_from')
		stop_index = elem.findtext('Seq-interval_to')
		gi = elem.findtext('Seq-interval_id/Seq-id/Seq-id_gi')
		return start_index, stop_index, gi
	
	def returnSeqLocPnt(self, elem):
		"""
		2008-07-29
		"""
		start_index = elem.findtext('Seq-point_point')
		stop_index = start_index
		gi = elem.findtext('Seq-point_id/Seq-id/Seq-id_gi')
		return start_index, stop_index, gi
	
	def return_location_list(self, elem):
		"""
		2008-07-29
			add one more case of getting location, "just one point"
			split codes into returnSeqLocInt() and returnSeqLocPnt()
		2008-07-28
			replace gi with gi_ls
			why +1 for all start and stop
		11-13-05 three versions
		"""
		start_ls = []
		stop_ls = []
		gi_ls = []
		seqint_elem_v1 = elem.find('Seq-loc/Seq-loc_mix/Seq-loc-mix')
		if seqint_elem_v1:
			for seq_loc_elem in seqint_elem_v1:	#multiple intervals or multi-points
				tmp_elem = seq_loc_elem.find('Seq-loc_int/Seq-interval')	#interval
				if tmp_elem:
					start_index, stop_index, gi = self.returnSeqLocInt(tmp_elem)
				else:
					tmp_elem = seq_loc_elem.find('Seq-loc_pnt/Seq-point')	#point
					if tmp_elem:
						start_index, stop_index, gi = self.returnSeqLocPnt(tmp_elem)
					else:
						start_index, stop_index, gi = -2, -2, None
				start_ls.append(int(start_index)+1)
				stop_ls.append(int(stop_index)+1)
				gi_ls.append(gi)
		else:
			tmp_elem = elem.find('Seq-loc/Seq-loc_int/Seq-interval')	#just one interval
			if tmp_elem:
				start_index, stop_index, gi = self.returnSeqLocInt(tmp_elem)
				start_ls.append(int(start_index)+1)
				stop_ls.append(int(stop_index)+1)
				gi_ls.append(gi)
			else:
				seqint_elem_v3 = elem.find('Seq-loc/Seq-loc_packed-int/Packed-seqint')	#a different multi-intervals
				if seqint_elem_v3:
					for tmp_elem in seqint_elem_v3:
						start_index, stop_index, gi = self.returnSeqLocInt(tmp_elem)
						start_ls.append(int(start_index)+1)
						stop_ls.append(int(stop_index)+1)
						gi_ls.append(gi)
				else:
					tmp_elem = elem.find('Seq-loc/Seq-loc_pnt/Seq-point')	#just one point
					if tmp_elem:
						start_index, stop_index, gi = self.returnSeqLocPnt(tmp_elem)
					else:	#run out of all choices
						start_index, stop_index, gi = -2, -2, None
					start_ls.append(int(start_index)+1)
					stop_ls.append(int(stop_index)+1)
					gi_ls.append(gi)
		return start_ls, stop_ls, gi_ls
	
	def returnGeneSegments(self, db, elem=None, gene_commentary=None, commentary_type=None):
		"""
		2012.5.15
			add argument commentary_type to stop replicating gene_commentary.gene_commentary_type
		2008-07-28
		"""
		start_ls, stop_ls, gi_ls = self.return_location_list(elem)
		gene_segments = []
		min_start = start_ls[0]
		max_stop = stop_ls[0]
		if commentary_type:
			gene_commentary_type = db.getGeneCommentaryType(commentary_type=commentary_type)
		else:
			gene_commentary_type = gene_commentary.gene_commentary_type
		for i in range(len(start_ls)):
			start = start_ls[i]
			stop = stop_ls[i]
			min_start_stop = min(start, stop)
			max_start_stop = max(start, stop)
			if min_start_stop < min_start:
				min_start = min_start_stop
			if max_start_stop > max_stop:
				max_stop = max_start_stop
			gi = gi_ls[i]
			gene_segment = GeneSegment(start=start, stop=stop, gi=gi, gene_commentary_type=gene_commentary_type)
			gene_segment.gene_commentary = gene_commentary
			gene_segments.append(gene_segment)
		passingdata = PassingData()
		passingdata.gene_segments = gene_segments
		passingdata.start = min_start
		passingdata.stop = max_stop
		return passingdata
	
	def returnGeneCommentary(self, elem, entrezgene_mapping, session, type=1):
		"""
		2012.5.15
			use self.db.getGeneCommentaryType() to get gene_commentary_type
		2008-07-29
			type is useless, GeneCommentary and GeneProduct merged
		2008-07-28
			type=1, GeneCommentary
			type=2, GeneProduct ( GeneCommentary in 'Gene-commentary_products')
		"""		
		commentary_type_id = int(elem.findtext('Gene-commentary_type'))
		commentary_type = elem.find('Gene-commentary_type').get('value')
		gene_commentary_type = self.db.getGeneCommentaryType(commentary_type=commentary_type)	#no passing of commentary_type_id
		
		if self.debug:
			sys.stderr.write("\t getting gene_commentary. type=%s.\n"%(gene_commentary_type.type))
		accession = elem.findtext('Gene-commentary_accession')
		version = elem.findtext('Gene-commentary_version')
		gi = elem.findtext('Gene-commentary_seqs/Seq-loc/Seq-loc_whole/Seq-id/Seq-id_gi')
		label = elem.findtext('Gene-commentary_label')
		text = elem.findtext('Gene-commentary_text')
		
		gene_commentary = GeneCommentary(label=label, accession=accession, version=version,\
										text=text, gi=gi)
		gene_commentary.gene_commentary_type = gene_commentary_type
		gene_commentary.gene_id = entrezgene_mapping.id
		#for the coordinates
		genomic_coords_elem = elem.find('Gene-commentary_genomic-coords')
		if genomic_coords_elem:
			if gene_commentary_type.type=='peptide':
				segment_commentary_type = 'CDS'
			elif gene_commentary_type.type=='mRNA':
				segment_commentary_type = 'exon'
			else:
				segment_commentary_type = None
			passingdata = self.returnGeneSegments(self.db, elem=genomic_coords_elem, gene_commentary=gene_commentary, \
												commentary_type=segment_commentary_type)
			if passingdata.gene_segments:
				gene_commentary.gene_segments = passingdata.gene_segments
				gene_commentary.start = passingdata.start
				gene_commentary.stop = passingdata.stop
		
		commentary_elems = elem.find('Gene-commentary_comment')
		if commentary_elems:
			comment_ls = []
			for commentary_elem in commentary_elems:
				commentary_type = commentary_elem.find('Gene-commentary_type').get('value')
				if self.debug:
					sys.stderr.write("\t\t getting Gene-commentary_comment type=%s.\n"%(commentary_type))
				if commentary_type=='comment':
					comment_summary_ls = []
					comment = commentary_elem.findtext('Gene-commentary_text')
					if comment:
						comment_summary_ls.append(comment)
					other_source_elem = commentary_elem.find('Gene-commentary_source')
					if other_source_elem:
						source_anchor = other_source_elem.findtext('Other-source/Other-source_anchor')
						if source_anchor:
							comment_summary_ls.append(source_anchor)
						source_url = other_source_elem.findtext('Other-source/Other-source_url')
						if source_url:
							comment_summary_ls.append(source_url)
					comment_ls.append(' '.join(comment_summary_ls))
			gene_commentary.comment = '|'.join(comment_ls)
		## to find further products (like peptide as products of mRNA)
		entrezgene_locus_products_elems = elem.find('Gene-commentary_products')
		if entrezgene_locus_products_elems:
			for entrezgene_locus_products_elem in entrezgene_locus_products_elems:
				if self.debug:
					sys.stderr.write("\t\t found gene_product commentary inside this gene_commentary.\n")
				sub_gene_commentary = self.returnGeneCommentary(entrezgene_locus_products_elem, entrezgene_mapping, session)
				#gene_product.gene_product = gene_commentary
				#session.add(gene_product)
				#session.flush()
				gene_commentary.gene_commentaries.append(sub_gene_commentary)
		return gene_commentary
		
	def find_peptide_info(self, elem, entrezgene_mapping):
		"""
		2008-07-28
			deprecated	
		11-13-05 for peptide
		"""
		entrezgene_mapping.cds_start, entrezgene_mapping.cds_stop, cds_loc_gi = self.return_location_list(elem)
		if cds_loc_gi:
			entrezgene_mapping.cds_gi = cds_loc_gi
		else:	#use the gi where the sequence is from
			entrezgene_mapping.cds_gi = cds_loc_gi
	
	def find_rna_info(self, elem, entrezgene_mapping):
		"""
		2008-07-28
			deprecated
		11-13-05 for rRNA and tRNA
		"""
		accession = elem.findtext('Gene-commentary_accession')
		label = elem.findtext('Gene-commentary_label')
		version = elem.findtext('Gene-commentary_version')
		if accession:
			entrezgene_mapping.mrna_acc_ver = accession + '.' + version
		elif label:
			entrezgene_mapping.mrna_acc_ver = label
		else:
			entrezgene_mapping.mrna_acc_ver = version
		mrna_gi = elem.findtext('Gene-commentary_seqs/Seq-loc/Seq-loc_whole/Seq-id/Seq-id_gi')
		entrezgene_mapping.mrna_start, entrezgene_mapping.mrna_stop, mrna_loc_gi = self.return_location_list(elem)
		if mrna_gi:
			entrezgene_mapping.mrna_gi = mrna_gi
		else:	#use the gi where the sequence is from
			entrezgene_mapping.mrna_gi = mrna_loc_gi
	
	def find_mrna_info(self, mrna_elem, entrezgene_mapping):
		"""
		2008-07-28
			deprecated
		11-13-05 for mRNA and its peptide
		"""
		#handle mRNA stuff
		entrezgene_mapping.mrna_acc_ver = mrna_elem.findtext('Gene-commentary_accession')\
			+'.'+mrna_elem.findtext('Gene-commentary_version')
		mrna_gi = mrna_elem.findtext('Gene-commentary_seqs/Seq-loc/Seq-loc_whole/Seq-id/Seq-id_gi')
		entrezgene_mapping.mrna_start, entrezgene_mapping.mrna_stop, mrna_loc_gi = self.return_location_list(mrna_elem)
		if mrna_gi:
			entrezgene_mapping.mrna_gi = mrna_gi
		else:	#use the gi where the sequence is from
			entrezgene_mapping.mrna_gi = mrna_loc_gi
		#handle the cds stuff
		cds_elem = mrna_elem.find('Gene-commentary_products/Gene-commentary')
		if cds_elem:
			entrezgene_mapping.cds_acc_ver = cds_elem.findtext('Gene-commentary_accession')\
				+'.'+cds_elem.findtext('Gene-commentary_version')
			cds_gi = cds_elem.findtext('Gene-commentary_seqs/Seq-loc/Seq-loc_whole/Seq-id/Seq-id_gi')
			entrezgene_mapping.cds_start, entrezgene_mapping.cds_stop, cds_loc_gi = self.return_location_list(cds_elem)
			if cds_gi:
				entrezgene_mapping.cds_gi = cds_gi
			else:	#use the gi where the sequence is from
				entrezgene_mapping.cds_gi = cds_loc_gi
		else:
			sys.stderr.write("\t gene_id: %s has no cds_elem.\n"%(entrezgene_mapping.gene_id))
	
	def return_datetime(self, date_elem):
		"""
		2010-12-15
			use returnAnyValueIfNothing() to deal with nothing returned.
		2008-07-28
		"""
		try:
			year = returnAnyValueIfNothing(date_elem.findtext('Date-std_year'), int)
			month= returnAnyValueIfNothing(date_elem.findtext('Date-std_month'), int)
			day = returnAnyValueIfNothing(date_elem.findtext('Date-std_day'), int)
			hour = returnAnyValueIfNothing(date_elem.findtext('Date-std_hour'))
			minute = returnAnyValueIfNothing(date_elem.findtext('Date-std_minute'))
			second = returnAnyValueIfNothing(date_elem.findtext('Date-std_second'))
			return datetime(year, month, day, hour, minute, second)
		except:
			traceback.print_exc()
			sys.stderr.write("Warning: %s.\n"%repr(sys.exc_info()))
			return None
	
	def parse_entrezgene_xml_file(self, session, filename):
		"""
		2010-12-14
			EntrezgeneMapping has been superceded by Gene.
		2008-07-28
		11-13-05 all xxx_from or xxx_to are computer indices starting from 0, so all +1
		"""
		counter = 0
		real_counter = 0
		for event, elem in ElementTree.iterparse(filename):
			if elem.tag == 'Entrezgene':
				gene_id = elem.findtext('Entrezgene_track-info/Gene-track/Gene-track_geneid')
				status_elem = elem.find('Entrezgene_track-info/Gene-track/Gene-track_status')
				status = status_elem.get('value')
				if status=='live':
					if self.debug:
						sys.stderr.write("gene_id=%s status: %s\n"%(gene_id, status))
					#entrezgene_mapping = entrezgene_mapping_attr()
					#entrezgene_mapping.gene_id = elem.findtext('Entrezgene_track-info/Gene-track/Gene-track_geneid')
					entrezgene_mapping = Gene.query.filter_by(ncbi_gene_id=gene_id).first()
					if entrezgene_mapping is None:
						entrezgene_mapping = Gene(ncbi_gene_id=gene_id)
						#sys.stderr.write("\t Warning: gene_id=%s not in table gene. EntrezgeneMapping entry skipped.\n"%(gi, gene_id))
						#continue
					entrezgene_mapping.date_created = self.return_datetime(elem.find('Entrezgene_track-info/Gene-track/Gene-track_create-date/Date/Date_std/Date-std'))
					entrezgene_mapping.date_updated = self.return_datetime(elem.find('Entrezgene_track-info/Gene-track/Gene-track_update-date/Date/Date_std/Date-std'))
					
					#2010-7-30 add 1000 to entrezgene_type_id because the original one could be 0 (unknown), not possible with id in mysql.
					entrezgene_type_id = int(elem.findtext('Entrezgene_type'))+1000
					entrezgene_type = EntrezgeneType.get(entrezgene_type_id)	#query the db to see if it exists or not
					if not entrezgene_type:
						entrezgene_type_value = elem.find('Entrezgene_type').get('value')
						entrezgene_type = EntrezgeneType(id=entrezgene_type_id, type=entrezgene_type_value)
						session.add(entrezgene_type)
						session.flush()
					entrezgene_mapping.entrezgene_type = entrezgene_type	#2010-7-30 not directly assigning the ID.
					
					entrezgene_mapping.tax_id = \
						elem.findtext('Entrezgene_source/BioSource/BioSource_org/Org-ref/Org-ref_db/Dbtag/Dbtag_tag/Object-id/Object-id_id')
					entrezgene_mapping.tax_id = int(entrezgene_mapping.tax_id)
					## see if there's chromosome information
					tmp_elem = elem.find('Entrezgene_source/BioSource/BioSource_subtype/SubSource')
					if tmp_elem!=None and tmp_elem.find('SubSource_subtype').get('value')=='chromosome':
						entrezgene_mapping.chromosome = tmp_elem.findtext('SubSource_name')
					else:
						tmp_elem = elem.find('Entrezgene_source/BioSource/BioSource_genome')
						if tmp_elem!=None:
							entrezgene_mapping.chromosome = tmp_elem.get('value')
					###  check location info
					el_elem = elem.find('Entrezgene_locus')
					if el_elem:	#some entries don't have Entrezgene_locus
						for el_sub_elem in el_elem:	#each el_sub_elem is Gene-commentary(omitted below) element.
							### 1st get genomic section
							el_sub_elem_type = el_sub_elem.find('Gene-commentary_type').get('value')
							if el_sub_elem.tag != 'Gene-commentary' or el_sub_elem_type!='genomic':
								continue
							seq_interval_elem = el_sub_elem.find('Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval')
							if not seq_interval_elem:	#no gi info, ignore, can't track genomic location
								if self.debug:
									sys.stderr.write("\t gene-id=%s has genomic commentary but no genomic position info.\n"%gene_id)
								continue
							gi = seq_interval_elem.findtext('Seq-interval_id/Seq-id/Seq-id_gi')
							gi =int(gi)
							if self.debug:
								sys.stderr.write("\t Entrezgene_locus genomic gi=%s.\n"%(gi))
							### ensure this gi is present in annot_assembly_table
							annot_assembly = AnnotAssembly.query.get(gi)
							if annot_assembly:
								#if self.is_gi_valid_in_annot_assembly_table(curs, gi, annot_assembly_table):
								if self.debug:
									sys.stderr.write("\t gi=%s found in annot_assembly_table.\n"%(gi))
								entrezgene_mapping.genomic_gi = gi
								entrezgene_mapping.genomic_accession = el_sub_elem.findtext('Gene-commentary_accession')
								entrezgene_mapping.genomic_version = int(el_sub_elem.findtext('Gene-commentary_version'))
								entrezgene_mapping.genomic_annot_assembly = annot_assembly	#2012.4.26 add the assembly
								strand_elem = seq_interval_elem.find('Seq-interval_strand/Na-strand')
								if strand_elem!=None:	#strand could be missing,
									#WATCH strand_elem has no children(len(strand_elem)=0), so 'if strand_elem:' doesn't work
									strand = strand_elem.get('value')
									if strand=="minus":
										entrezgene_mapping.strand = '-1'
									elif strand=="plus":
										entrezgene_mapping.strand = '+1'
									else:
										entrezgene_mapping.strand = strand
									if self.debug:
										sys.stderr.write("\t entrezgene_mapping.strand=%s.\n"%entrezgene_mapping.strand)
								entrezgene_mapping.start = int(seq_interval_elem.findtext('Seq-interval_from'))+1
								entrezgene_mapping.stop = int(seq_interval_elem.findtext('Seq-interval_to'))+1
								session.add(entrezgene_mapping)
								session.flush()
								### 2nd check the products of the genomic region
								entrezgene_locus_products_elems = el_sub_elem.find('Gene-commentary_products')
								if entrezgene_locus_products_elems:
									for entrezgene_locus_products_elem in entrezgene_locus_products_elems:
										gene_commentary_type = entrezgene_locus_products_elem.find('Gene-commentary_type').get('value')
										if self.debug:
											sys.stderr.write("\t entrezgene_locus_products_elem type: %s\n"%gene_commentary_type)
										gene_commentary = self.returnGeneCommentary(entrezgene_locus_products_elem, entrezgene_mapping, session)
										session.add(gene_commentary)
										session.flush()
								## some entrezgene_locus have commentary_comment
								entrezgene_locus_genomic_commentary_comment_elems = el_sub_elem.find('Gene-commentary_comment')
								if entrezgene_locus_genomic_commentary_comment_elems:
									for tmp_elem in entrezgene_locus_genomic_commentary_comment_elems:
										gene_commentary = self.returnGeneCommentary(tmp_elem, entrezgene_mapping, session)
										session.add(gene_commentary)
										session.flush()
								
								#submit to database
								#self.submit_to_entrezgene_mapping_table(curs, table, entrezgene_mapping)
								real_counter += 1
							else:
								sys.stderr.write("\t Warning: genomic_gi=%s, not in annot_assembly_table gene_id=%s\n"%(gi, gene_id))
				#release memory
				elem.clear()
				counter += 1
			if counter%5000==0:
				session.flush()
				#session.clear()	#2010-6-22 gone in 0.6
				session.expunge_all()	# 2010-6-22 sqlalchemy 0.6
			if self.report and counter%2000==0:
				sys.stderr.write("%s\t%s/%s"%('\x08'*20, counter, real_counter))
		if self.report:
			sys.stderr.write("%s\t%s/%s\n"%('\x08'*20, counter, real_counter))
	
	def run(self):
		"""
		11-13-05 
			--db_connect()
			--parse_entrezgene_xml_file()
				--is_gi_valid_in_annot_assembly_table()
				--find_info_dict()
					--return_location_list()
				--submit_to_entrezgene_mapping_table()
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.inputfiles))
		db = GenomeDatabase(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)	#2010-6-22
		self.db = db
		session = db.session
		if not self.debug:	#in debug mode, no transaction, auto-commit
			session.begin()
		
		#conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		for f in self.inputfiles:
			sys.stderr.write("%d/%d:\t%s\n"%(self.inputfiles.index(f)+1,len(self.inputfiles),f))
			self.parse_entrezgene_xml_file(session, f)
		if not self.debug:
			if self.commit:
				session.commit()
				#curs.execute("end")
			else:
				session.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = GeneASNXML2gene_mapping
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:a:cbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'sequence'
	table = 'entrezgene_mapping'
	annot_assembly_table = 'annot_assembly'
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
		elif opt in ("-t"):
			table = arg
		elif opt in ("-a"):
			annot_assembly_table = arg
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if len(args)>0:
		instance = GeneASNXML2gene_mapping(hostname, dbname, schema, args, table, annot_assembly_table, \
			commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
	"""
