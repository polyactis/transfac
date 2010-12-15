#!/usr/bin/env python
"""
Examples:
	#populate gene database in postgresql
	UpdateGeneDB.py -u crocea -c
	
	# 2010-12-14 add gene info into schema vervetdb.genome. files stored in /usr/local/research_data/NCBI/
	UpdateGeneDB.py -d vervetdb -k genome -u yh -c -v postgresql -t /usr/local/research_data/NCBI/
	
	#populate gene database in mysql
	UpdateGeneDB.py -v mysql -u yh -z papaya -d genome -k "" -c
	
Description:
	Program to fill the genome db with Gene infomation files, which are downloaded from the web.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt, gzip
import warnings, traceback, gc, subprocess
from UpdateGenomeDB import UpdateGenomeDB
from pymodule.GenomeDB import GenomeDatabase, Gene, Gene2go

class UpdateGeneDB(UpdateGenomeDB):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgres', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['graphdb', 'd', 1, 'database name', ],\
							('schema', 1, ): ['genome', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('tmp_dir',1,):['/tmp/', 't', 1, 'Temporary directory to hold files downloaded from web'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-06
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	gene_info_ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz'
	
	gene2go_ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'
	
	def submitGeneInfo(self, gene_info_fname, session, input_fname_type=1):
		"""
		2008-07-08
			input_fname_type=1: gene_info_fname=gene_info.gz
			input_fname_type=2: gene_info_fname=gene2go.gz
		"""
		sys.stderr.write("Submitting %s into db ... \n"%gene_info_fname)
		if gene_info_fname[-2:]=='gz':
			inf = gzip.open(gene_info_fname, 'r')
		else:
			inf = open(gene_info_fname, 'r')
		reader = csv.reader(inf, delimiter='\t')
		reader.next()	#skip the header
		counter = 0
		real_counter = 0
		for row in reader:
			for i in range(len(row)):
				if row[i]=='-':
					row[i] = None
			counter += 1
			if input_fname_type==1:
				tax_id, ncbi_gene_id, gene_symbol, locustag, synonyms, dbxrefs, chromosome, map_location, description, \
				type_of_gene, symbol_from_nomenclature_authority, full_name_from_nomenclature_authority, nomenclature_status,\
				other_designations, modification_date = row
				
				gene = Gene(tax_id=tax_id, ncbi_gene_id=ncbi_gene_id, gene_symbol=gene_symbol, locustag=locustag, \
						synonyms=synonyms, dbxrefs=dbxrefs, chromosome=chromosome, map_location=map_location, \
						description=description, \
						type_of_gene=type_of_gene, symbol_from_nomenclature_authority=symbol_from_nomenclature_authority, \
						full_name_from_nomenclature_authority=full_name_from_nomenclature_authority, nomenclature_status=nomenclature_status,\
						other_designations=other_designations, modification_date=modification_date)
				session.add(gene)
				real_counter += 1
			elif input_fname_type==2:
				tax_id, ncbi_gene_id, go_id, evidence, go_qualifier, go_description, pubmed_ids, category = row[:8]
				gene = Gene.query.filter_by(ncbi_gene_id=ncbi_gene_id).first()
				if gene is not None:
					gene2go = Gene2go(tax_id=tax_id, go_id=go_id, evidence=evidence, go_qualifier=go_qualifier, \
								go_description=go_description, pubmed_ids=pubmed_ids, category=category)
					gene2go.gene = gene
					session.add(gene2go)
					real_counter += 1
			if counter%5000==0:
				session.flush()
				#session.clear()	#2010-12-14 gone in 0.6
				session.expunge_all()	# 2010-12-14 sqlalchemy 0.6
				sys.stderr.write("%s\t%s\t%s(skipped)"%('\x08'*40, counter, counter-real_counter))
		del reader, inf
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2010-12-14
			add "db.setup(create_tables=False)" to setup db structure
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db = GenomeDatabase(drivername=self.drivername, username=self.db_user,
				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		#session.begin()
		gene_info_fname, gene2go_fname = self.getInputFileList([self.gene_info_ftp_url, self.gene2go_ftp_url], self.tmp_dir)
		#gene_info_fname = '/usr/local/research_data/ncbi/gene_2008_07_06/DATA/gene_info.gz'
		self.submitGeneInfo(gene_info_fname, session)
		#gene2go_fname = '/usr/local/research_data/ncbi/gene_2008_07_06/DATA/gene2go.gz'
		self.submitGeneInfo(gene2go_fname, session, input_fname_type=2)
		if self.commit:
			session.flush()
			session.commit()
			#db.transaction.commit()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()
			#db.transaction.rollback()
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = UpdateGeneDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
