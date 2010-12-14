#!/usr/bin/env python
"""
Examples:
	#postgres
	UpdateGenomeDB.py -g at -c
	
	#mysql
	UpdateGenomeDB.py -v mysql -u yh -z papaya -d genome -k "" -g at -c
	
Description:
	A wrapper program to fill database with genome sequences. Call chromosome_fasta2db for each file downloaded.

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
							('organism', 0, ): [None, 'g', 1, '2-letter abbreviation for organism. Optional, if specified, only sequence from this organism would be extracted.', ],\
							('tmp_dir',1,):['/tmp/', 't', 1, 'Temporary directory to hold files downloaded from web'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	short_org2chr_file_ls = {'at':['ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_I/NC_003070.fna',\
								'ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_II/NC_003071.fna',\
								'ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_III/NC_003074.fna',\
								'ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_IV/NC_003075.fna',\
								'ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/CHR_V/NC_003076.fna']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-06
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
	def getInputFileList(self, file_source_ls, tmp_dir):
		"""
		2008-07-06
		"""
		sys.stderr.write("Getting input file list ...")
		command = 'wget'
		file_list = []
		for file_source in file_source_ls:
			filename = os.path.basename(file_source)
			local_filename = os.path.join(tmp_dir, filename)
			if os.path.isfile(local_filename):
				sys.stderr("%s exists. Skip.\n"%local_filename)
			else:
				cmd_p = subprocess.Popen([command, '-P', tmp_dir, file_source], stderr=sys.stderr, stdout=sys.stdout)
				cmd_p.wait()
			file_list.append(local_filename)
		sys.stderr.write("Done.\n")
		return file_list
	
	def run(self):
		file_list = self.getInputFileList(self.short_org2chr_file_ls[self.organism], self.tmp_dir)
		for filename in file_list:
			cf2db_ins = chromosome_fasta2db(drivername=self.drivername, hostname=self.hostname, dbname=self.dbname, schema=self.schema, \
										db_user=self.db_user, db_passwd=self.db_passwd, inputfiles=[filename], \
										organism=self.organism, commit=self.commit, report=self.report, debug=self.debug)
			cf2db_ins.run()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = UpdateGenomeDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()