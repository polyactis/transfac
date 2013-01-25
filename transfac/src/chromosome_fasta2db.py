#!/usr/bin/env python
"""

Examples:
	chromosome_fasta2db.py -c -i *.fa.gz
	
	#for mysql
	chromosome_fasta2db.py -v mysql -u yh -z papaya -d genome -k "" -c -i *.fa.gz
	chromosome_fasta2db.py -v mysql -u yh -z papaya -d genome -k "" -c -i /usr/local/research_data/ncbi/genomes_2008_07_29/mitochondrion1.genomic.fna.gz -b -r
	
	#take only A. thaliana's mitochondrion
	chromosome_fasta2db.py -v mysql -u yh -z papaya -d genome -k "" -c -i mitochondrion1.genomic.fna.gz -b -r -g at
	
	# 2011-7-7 put BACs into db
	chromosome_fasta2db.py -u yh -k genome -d vervetdb -v postgresql -c -i script/vervet/data/ref/BAC/BAC.accession.fasta
		-y 3 -b -r -s BAC -g "Chlorocebus aethiops"
	
	# 2011-7-7 put 1000 scaffolds (order in the file) into db
	chromosome_fasta2db.py -u yh -k genome -d vervetdb -v postgresql -c 
	-i script/vervet/data/Draft_June_2011/supercontigs/supercontigs.fasta -y 2  -r -s Scaffold -g "Chlorocebus aethiops" -x 1000
	
Description:
	Parse the chromosome sequence files (fasta format) downloaded from
		ftp://ftp.ncbi.nih.gov/genomes/ORGANISM/Assembled_chromosomes/
		ftp://ftp.ncbi.nih.gov/refseq/release/
	INPUT_FILEs are *.fa.gz (fasta format). The program detects whether it's gzipped or not
		based on the filename (like *gz or not). It could contain multiple fasta blocks.
		Taxonomy ID of the organism is detected thru the fasta header by FigureOutTaxID
			based on data stored in tables of schema taxonomy in the same postgres db, on the same host.
	
	Results are stored in 2 tables in the database, annot_assembly and raw_sequence.
		The chromosome sequence is snipped into 10kb chunks and dumped into raw_sequence_table.
		If annot_assembly has the same entry based on the unique constraint, this entry would be ignored. 
	
"""

import sys, getopt, os, re, gzip
sys.path += [os.path.join(os.path.expanduser('~/script'))]
from annot.bin.codense.common import db_connect, org_short2long, org2tax_id
from pymodule.GenomeDB import *
from pymodule import PassingData
from pymodule.utils import FigureOutTaxID
#from pymodule.db import db_connect

class annot_assembly_attr:
	"""
	2008-07-27
		deprecated
	"""
	def __init__(self):
		self.gi = None
		self.acc_ver = None
		self.acc = None
		self.version = None
		self.tax_id = None
		self.chromosome = None
		self.start = None
		self.stop = None
		self.orientation = None
		self.sequence = None
		self.raw_sequence_start_id = None
		self.seq_type = None
		self.comment = None

class chromosome_fasta2db:
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgres', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['graphdb', 'd', 1, 'database name', ],\
							('schema', 1, ): ['genome', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('inputfiles', 1, ):[None, 'i', 1, 'comma-separated input filenames, or a list of files'],\
							('organism', 0, ): [None, 'g', 1, '2-letter abbreviation for organism. Optional, if specified, only sequence from this organism would be extracted.'],\
							('sequence_type_name', 1, ):['chromosome sequence', 's', 1, 'short name in table sequence_type'],\
							('run_type', 1, int):[1, 'y', 1, 'run type. 1: genBank fasta files. 2: scaffolds from WUSTL. 3: fully sequenced vervet BACs. 4: fully-sequenced vervet BAC ends'],\
							('maxNoOfFastaRecords', 1, int):[500, 'x', 1, 'maximum number of fasta records to be inserted (in whatever order)'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-07-27
			use option_default_dict
		2008-07-06
			use the firstline (header) of the fasta file to extract which chromosome. using filename is unreliable.
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		if type(self.inputfiles)==str:
			self.inputfiles = self.inputfiles.split(',')
		
		self.FigureOutTaxID_ins = FigureOutTaxID(db_user=self.db_user,
								db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname)
		if self.organism is not None:
			if org_short2long(self.organism):
				self.tax_id = org2tax_id(org_short2long(self.organism))
			else:
				self.tax_id = self.FigureOutTaxID_ins.returnTaxIDGivenSentence(self.organism)
		
		else:
			self.tax_id=None
		
		#self.p_chromosome = re.compile(r'[a-zA-Z]+_chr(\w+).fa')
		self.p_chromosome = re.compile(r'chromosome (\w+)[,\n\r]?')	#the last ? means [,\n\r] is optional
		self.p_acc_ver = re.compile(r'(\w+)\.(\d+)')
		
		self.parseFastaDescriptionDict = {1: self.parseFastaDescriptionForGenBank, 2: self.parseFastaDescriptionForWUSTLVervetScaffolds,
										3: self.parseFastaDescriptionForFullVervetBACs}
	
	def saveRawSequence(self, session, seq_to_db, passingdata, aa_attr_instance):
		"""
		2010-12-17
			RawSequence.annot_assembly is a foreign key element now.
		2008-07-29
			to store one sequence segment
		"""
		passingdata.current_stop = passingdata.current_start+len(seq_to_db)-1
		raw_sequence = RawSequence(start=passingdata.current_start, stop=passingdata.current_stop, sequence=seq_to_db)
		raw_sequence.annot_assembly = aa_attr_instance
		session.add(raw_sequence)
		if not passingdata.raw_sequence_initiated:
			session.flush()	# 2010-12-17 to get raw_sequence.id
			passingdata.raw_sequence_initiated = True
			aa_attr_instance.raw_sequence_start_id = raw_sequence.id
		passingdata.current_start += len(seq_to_db)
	
	def parseFastaDescriptionForGenBank(self, descriptionLine=None, FigureOutTaxID_ins=None):
		"""
		2011-7-6
			
		"""
		"""
		possible header lines:
		
		>gi|51511461|ref|NC_000001.8|NC_000001 Homo sapiens chromosome 1, complete sequence
		>gi|186497660|ref|NC_003070.6| Arabidopsis thaliana chromosome 1, complete sequence
		>gi|26556996|ref|NC_001284.2| Arabidopsis thaliana mitochondrion, complete genome
		>gi|115442598|ref|NC_008394.1| Oryza sativa (japonica cultivar-group) genomic DNA, chromosome 1
		"""
		header = descriptionLine[1:-1]	#discard '>' and '\n'
		header = header.split('|')
		_tax_id = FigureOutTaxID_ins.returnTaxIDGivenSentence(header[4])
		
		if self.p_chromosome.search(header[4]) is not None:
			chromosome = self.p_chromosome.search(header[4]).groups()[0]
		elif header[4].find('mitochondrion')!=-1:
			chromosome = 'mitochondrion'
		elif header[4].find('chloroplast')!=-1:
			chromosome = 'chloroplast'
		else:	#something else, take the whole before ','
			chromosome = header[4].split(',')[0]
		gi = int(header[1])
		acc_ver = header[3]
		comment = header[4]
		return PassingData(tax_id=_tax_id, gi=gi, comment=comment, acc_ver=acc_ver, chromosome=chromosome)

	def parseFastaDescriptionForWUSTLVervetScaffolds(self, descriptionLine=None, FigureOutTaxID_ins=None):
		"""
		2011-7-6
			
		"""
		"""
		possible header lines:
		>Contig0  12652774 13406928

		"""
		header = descriptionLine[1:-1]	#discard '>' and '\n'
		header = header.split()
		chromosome = header[0]	#contig name is taken as chromosome
		"""
		p_chromosome = re.compile(r'Contig(\d+)')
		if p_chromosome.search(header[0]) is not None:
			chromosome = p_chromosome.search(header[0]).groups()[0]
		else:
			chromosome = None
		"""
		gi = None
		acc_ver = None
		comment = None
		return PassingData(tax_id=None, gi=gi, comment=comment, acc_ver=acc_ver, chromosome=chromosome)
	
	
	def parseFastaDescriptionForFullVervetBACs(self, descriptionLine=None, FigureOutTaxID_ins=None):
		"""
		2011-7-6
			
		possible header lines:
			
		>gi|285026568|gb|AC239257.2| Chlorocebus aethiops chromosome UNK clone CH252-270J24, WORKING DRAFT SEQUENCE, 2 unordered pieces
		>gi|281332227|gb|AC238852.3| Chlorocebus aethiops BAC clone CH252-133A18 from chromosome 3, complete sequence
		>gi|285002488|gb|AC239185.3| Chlorocebus aethiops BAC clone CH252-404N12 from chromosome unknown, complete sequence
		
		
		"""
		header = descriptionLine[1:-1]	#discard '>' and '\n'
		header = header.split('|')
		_tax_id = None
		p_chromosome = re.compile(r'UNK clone ([^,]+),')	# 1st type of clone description
		p2_chromosome = re.compile(r'clone ([^,]+),')	# 2nd type of clone description
		
		if p_chromosome.search(header[4]) is not None:
			chromosome = p_chromosome.search(header[4]).groups()[0]
		else:
			if p2_chromosome.search(header[4]) is not None:
				chromosome = p2_chromosome.search(header[4]).groups()[0]
			else:
				chromosome = None
		gi = int(header[1])
		acc_ver = header[3]
		comment = header[4]
		return PassingData(tax_id=_tax_id, gi=gi, comment=comment, acc_ver=acc_ver, chromosome=chromosome)
	
	def parse_chromosome_fasta_file(self, session, filename, gzipped, tax_id=None, chunk_size=10000, \
								sequence_type_name='chromosome sequence', run_type=1,
								maxNoOfFastaRecords=500):
		"""
		2011-7-10
			add argument maxNoOfFastaRecords: the max number of fasta records before quitting
		2011-7-6
			add argument run_type
				1: chromosome sequences from NCBI genbank
				2: vervet scaffolds from WUSTL
				3: full vervet BACs from McGill
		2010-12-15
			fix a bug that _tax_id shall be used in query AnnotAssembly.
			This bug caused the db redundancy check to fail.
		2010-12-15
			if entry already exists in AnnotAssembly, skip it.
		2008-07-29
			figure out tax_id via FigureOutTaxID
			filename could contain multiple fasta blocks
		2008-07-27
			change to use data structures from GenomeDB.py
		2008-07-06
			use the firstline (header) of the fasta file to extract which chromosome. using filename is unreliable.
		"""
		if gzipped:
			inf = gzip.open(filename, 'r')
		else:
			inf = open(filename, 'r')
		
		line = inf.readline()
		new_fasta_block = 1	#'line' is not enough to stop the 'while' loop. after the file reading is exhausted by "for line in inf:", 'line' still contains the stuff from the last line.
		no_of_fasta_blocks = 0
		while line and new_fasta_block:
			new_fasta_block = 0	#set it to 0, assuming only one fasta block, change upon new fasta block
			if line[0]!='>':	#not fasta block header
				for line in inf:	#exhaust this fasta block as it's not what's wanted.
					if line[0]=='>':
						new_fasta_block = 1
						break	#start from while again
				continue
			headerData = self.parseFastaDescriptionDict[run_type](line, self.FigureOutTaxID_ins)
			
			if tax_id is not None and tax_id!=headerData.tax_id:
				sys.stderr.write("tax_id (%s) not matching the one given (%s). Ignore.\n"%(headerData.tax_id, tax_id))
				line = inf.readline()
				new_fasta_block = 1
				continue
			
			chromosome = headerData.chromosome
			sequence_type = SequenceType.query.filter_by(type=sequence_type_name).first()
			start = 1
			aa_attr_instance = AnnotAssembly.query.filter_by(chromosome=chromosome).filter_by(tax_id=tax_id).filter_by(start=start).\
				filter_by(sequence_type_id=sequence_type.id).first()
			if aa_attr_instance and aa_attr_instance.raw_sequence_start_id is not None:
				# if raw sequences have been associated with this AnnotAssembly and 
				sys.stderr.write("raw sequences have been associated with this AnnotAssembly (tax_id %s, chr=%s, start=%s). Ignore.\n"%\
								(tax_id, chromosome, start))
				line = inf.readline()
				new_fasta_block = 1
				continue
			if aa_attr_instance is None:
				aa_attr_instance = AnnotAssembly()
				aa_attr_instance.gi = headerData.gi
				aa_attr_instance.acc_ver = headerData.acc_ver
				if aa_attr_instance.acc_ver and self.p_acc_ver.search(aa_attr_instance.acc_ver):
					aa_attr_instance.accession, aa_attr_instance.version = self.p_acc_ver.search(aa_attr_instance.acc_ver).groups()
					aa_attr_instance.version = int(aa_attr_instance.version)
				else:
					aa_attr_instance.accession = None
					aa_attr_instance.version = None
				aa_attr_instance.tax_id = tax_id
				if self.debug:
					sys.stderr.write("tax_id=%s for %s.\n"%(aa_attr_instance.tax_id, line))
				
				
				aa_attr_instance.chromosome = chromosome
				aa_attr_instance.start = 1
				#aa_attr_instance.raw_sequence_start_id = self.get_current_max_raw_sequence_id(curs, raw_sequence_table)+1
				aa_attr_instance.sequence_type_id = sequence_type.id
				aa_attr_instance.comment = headerData.comment	#store this whole thing for future reference
			
			passingdata = PassingData()
			passingdata.current_start = 1
			passingdata.raw_sequence_initiated = False
			seq = ''
			for line in inf:
				if line[0]=='>':
					if seq:	#last segment from the previous fasta block
						self.saveRawSequence(session, seq, passingdata, aa_attr_instance)
						seq = ''	#set to nothing to avoid saving one more RawSequence
					new_fasta_block = 1
					break	#start from while again
				
				seq += line.strip()
				if len(seq)>=chunk_size:
					seq_to_db = seq[:chunk_size]
					self.saveRawSequence(session, seq_to_db, passingdata, aa_attr_instance)
					seq = seq[chunk_size:]	#remove the one already in db
					if self.report:
						sys.stderr.write("%s\t%s\t%s"%('\x08'*20, no_of_fasta_blocks, passingdata.current_start/chunk_size+1))
			if seq:	# last segment from last line
				self.saveRawSequence(session, seq, passingdata, aa_attr_instance)
			aa_attr_instance.stop = passingdata.current_stop
			session.add(aa_attr_instance)
			session.flush()
			no_of_fasta_blocks += 1
			if no_of_fasta_blocks>=maxNoOfFastaRecords:
				break
		if self.report:
			sys.stderr.write("Number of fasta blocks: %s.\n"%(no_of_fasta_blocks))
		del inf
	
	
	def run(self):
		"""
		2008-07-27
			
			--GenomeDatabase
			--parse_chromosome_fasta_file()
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.inputfiles))
		db = GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		session.begin()
		for f in self.inputfiles:
			sys.stderr.write("%d/%d:\t%s\n"%(self.inputfiles.index(f)+1,len(self.inputfiles),f))
			if f[-2:]=='gz':	#turn on the gzip flag based on the filename extension
				gzipped=1
			else:
				gzipped=0
			self.parse_chromosome_fasta_file(session, f, gzipped, self.tax_id, sequence_type_name=self.sequence_type_name, \
									run_type=self.run_type, maxNoOfFastaRecords=self.maxNoOfFastaRecords)
		if self.commit:
			session.commit()
		else:
			session.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = chromosome_fasta2db
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:a:g:pcr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'sequence'
	annot_assembly_table = 'annot_assembly'
	raw_sequence_table = 'raw_sequence'
	organism = 'hs'
	gzipped = 0
	commit = 0
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
			annot_assembly_table = arg
		elif opt in ("-a"):
			raw_sequence_table = arg
		elif opt in ("-g"):
			organism = arg
		elif opt in ("-p"):
			gzipped = 1
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-r"):
			report = 1
	if len(args)>0:
		instance = chromosome_fasta2db(hostname, dbname, schema, args, annot_assembly_table, \
			raw_sequence_table, organism, gzipped, commit, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
	"""