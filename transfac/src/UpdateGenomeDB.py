#!/usr/bin/env python
"""
Examples:
	#postgres
	UpdateGenomeDB.py -g at -c
	
	# 2010-12-14 add human genome sequences into vervetdb.genome.
	UpdateGenomeDB.py -g hs -d vervetdb -k genome -u yh -c -v postgresql -t /usr/local/research_data/NCBI/hs
	
	# 2010-12-14 add chimp genome sequences into vervetdb.genome.
	UpdateGenomeDB.py -g ptr -d vervetdb -k genome -u yh -c -v postgresql -t /usr/local/research_data/NCBI/ptr
	
	#mysql
	UpdateGenomeDB.py -v mysql -u yh -z papaya -d genome -k "" -g at -c
	
Description:
	A wrapper program to fill database with genome sequences. Call chromosome_fasta2db for each file downloaded.
	
	The sequence file from NCBI contains the scientific name of the species,
		which would be used to find out the tax-id through querying tables in the taxonomy schema.
	The organism argument of the program is only used to download sequences from web.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback, gc, subprocess
from chromosome_fasta2db import chromosome_fasta2db

class UpdateGenomeDB(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgres', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['graphdb', 'd', 1, 'database name', ],\
							('schema', 1, ): ['genome', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('organism', 0, ): [None, 'g', 1, '2-letter abbreviation for organism. used to download sequence files.', ],\
							('tmp_dir',1,):['/tmp/', 't', 1, 'Temporary directory to hold files from web. Will be created if not present'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	short_org2chr_file_ls = {'at':['ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_I/NC_003070.fna',\
								'ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_II/NC_003071.fna',\
								'ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_III/NC_003074.fna',\
								'ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_IV/NC_003075.fna',\
								'ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_V/NC_003076.fna',\
								'ftp://ftp.ncbi.nih.gov/refseq/release/mitochondrion/mitochondrion1.genomic.fna.gz',\
								'ftp://ftp.ncbi.nih.gov/refseq/release/plastid/plastid1.genomic.fna.gz']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-06
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	@classmethod
	def getInputFileList(cls, file_source_ls, tmp_dir):
		"""
		2011-1-22
			become classmethod
		2010-12-14
			create new directory if tmp_dir is not present
		2008-07-06
		"""
		sys.stderr.write("Getting input file list ...")
		command = 'wget'
		file_list = []
		for file_source in file_source_ls:
			filename = os.path.basename(file_source)
			local_filename = os.path.join(tmp_dir, filename)
			if os.path.isfile(local_filename):
				sys.stderr.write("%s exists. Skip.\n"%local_filename)
			else:
				if not os.path.isdir(tmp_dir):
					os.makedirs(tmp_dir)
				cmd_p = subprocess.Popen([command, '-P', tmp_dir, file_source], stderr=sys.stderr, stdout=sys.stdout)
				cmd_p.wait()
			file_list.append(local_filename)
		sys.stderr.write("Done.\n")
		return file_list
	
	def run(self):
		"""
		2010-12-15
			add URLs for human, chimp, macaque genome sequences.
		"""
		#2010-12-14 URl for the human genome sequences
		speciesID = 'hs'
		self.short_org2chr_file_ls[speciesID] = []
		for chr in range(1,23)+['X', 'Y','MT']:
			self.short_org2chr_file_ls[speciesID].\
				append('ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh37.p2_chr%s.fa.gz'%chr)
			# 2010-12-15 the celera version is not needed
			#self.short_org2chr_file_ls[speciesID].\
			#	append('ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/H_sapiens/Assembled_chromosomes/seq/hs_alt_Hs_Celera_chr%s.fa.gz'%chr)
		
		#2010-12-14 URl for the rhesus macaque genome sequences
		speciesID = 'mm'
		self.short_org2chr_file_ls[speciesID] = []
		for chr in range(1,21)+['X',]:
			self.short_org2chr_file_ls[speciesID].\
				append('ftp://ftp.ncbi.nlm.nih.gov/genomes/Macaca_mulatta/Assembled_chromosomes/seq/mmu_ref_Mmul_051212_chr%s.fa.gz'%chr)
		self.short_org2chr_file_ls[speciesID].\
				append('ftp://ftp.ncbi.nlm.nih.gov/genomes/Macaca_mulatta/Assembled_chromosomes/seq/mmu_chrMT.fa.gz')
		
		#2010-12-14 URl for the chimp
		speciesID = 'ptr'
		self.short_org2chr_file_ls[speciesID] = []
		for chr in [1]+range(3,23)+['X','Y', '2A', '2B']:	#chimp doesn't have chr 2.
			self.short_org2chr_file_ls[speciesID].\
				append('ftp://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/ptr_ref_chr%s.fa.gz'%chr)
		# 2010-12-15 an alternative chimp Y chromosome assembly (cause confusing in terms of gene coordinates.)
		#self.short_org2chr_file_ls[speciesID].\
		#		append('ftp://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/ptr_alt_CCYSCv1_chrY.fa.gz')
		
		#2010-12-17 URl for the orangutan
		speciesID = 'pa'
		self.short_org2chr_file_ls[speciesID] = []
		for chr in range(1,23)+['X','Y', '2A', '2B']:	# no chr2
			self.short_org2chr_file_ls[speciesID].\
				append("ftp://ftp.ncbi.nlm.nih.gov/genomes/Pongo_abelii/Assembled_chromosomes/seq/pab_ref_P_pygmaeus_2.0.2_chr%s.fa.gz"%chr)
		self.short_org2chr_file_ls[speciesID].\
			append("ftp://ftp.ncbi.nlm.nih.gov/genomes/Pongo_abelii/Assembled_chromosomes/seq/pab_chrMT.fa.gz")
		
		# 2010-12-17 "unlocalized" or "unplaced" would mess up the chromosome inference code. which will regard them as chromosome 1.
		#self.short_org2chr_file_ls[speciesID].\
		#	append("ftp://ftp.ncbi.nlm.nih.gov/genomes/Pongo_abelii/Assembled_chromosomes/seq/pab_ref_P_pygmaeus_2.0.2_unlocalized.fa.gz")
		#self.short_org2chr_file_ls[speciesID].\
		#	append("ftp://ftp.ncbi.nlm.nih.gov/genomes/Pongo_abelii/Assembled_chromosomes/seq/pab_ref_P_pygmaeus_2.0.2_unplaced.fa.gz")
		
		# 2010-12-17 URL for marmoset
		speciesID = 'cj'
		self.short_org2chr_file_ls[speciesID] = []
		for chr in range(1,23)+['X','Y', ]:
			self.short_org2chr_file_ls[speciesID].\
				append("ftp://ftp.ncbi.nlm.nih.gov/genomes/Callithrix_jacchus/Assembled_chromosomes/seq/cja_ref_Callithrix_jacchus-3.2_chr%s.fa.gz"%chr)
		
		file_list = self.getInputFileList(self.short_org2chr_file_ls[self.organism], self.tmp_dir)
		for filename in file_list:
			if os.path.isfile(filename):
				#2010-12-14 stop passing self.organism to chromosome_fasta2db because the latter can parse tax_id out from files.
				# and two-letter organism symbols start to be non-unique.
				cf2db_ins = chromosome_fasta2db(drivername=self.drivername, hostname=self.hostname, dbname=self.dbname, schema=self.schema, \
											db_user=self.db_user, db_passwd=self.db_passwd, inputfiles=[filename], \
											organism=None, commit=self.commit, report=self.report, debug=self.debug)
				cf2db_ins.run()
			else:
				sys.stderr.write("Warning: file %s doesn't exist. Skipped.\n"%(filename))

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = UpdateGenomeDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
