#!/usr/bin/env python
"""
2008-07-08
A wrapper for the genome database built on top of sqlalchemy. derived from variation.src.db
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import sqlalchemy, threading
from sqlalchemy.engine.url import URL
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey
from sqlalchemy.orm.session import Session
from sqlalchemy.orm import mapper, relation

from pymodule.db import Database, TableClass


class Gene(TableClass):
	pass

class Gene2go(TableClass):
	pass

class Gene_symbol2id(TableClass):
	pass

class README(TableClass):
	pass

class GenomeDatabase(Database):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgres', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['localhost', 'z', 1, 'hostname of the db server', ],\
							('database', 1, ):['graphdb', 'd', 1, '',],\
							('schema', 1, ): ['genome', 'k', 1, 'database schema name', ],\
							('username', 1, ):['crocea', 'u', 1, 'database username',],\
							('password',0, ):[None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, 'o', 1, 'database port number'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self._threadlocal = threading.local()
		self.tables = {}
		self.mappers = {}
		self._engine = None
	
	def _setup_tables(self, metadata, tables):
		"""
		2008-07-09
			specify the schema (postgres) in Table()
		Map the database structure to SQLAlchemy Table objects
		"""
		table_ls = ['gene', 'gene2go', 'gene_symbol2id', 'readme']
		for table_name in table_ls:
			tables[table_name] = Table(table_name, metadata, autoload=True, schema=getattr(self, 'schema', None))	#2008-07-09 specify the schema (postgres)
	
	#_setup_tables = classmethod(_setup_tables)
	
	def _setup_mappers(cls, tables, mappers):
		"""
		2008-07-08
		"""
		standalone_table_tuple_ls = [('gene', Gene), ('gene2go', Gene2go), ('gene_symbol2id', Gene_symbol2id), ('readme', README)]
		for table_name, table_class in standalone_table_tuple_ls:
			mappers[table_name] = mapper(table_class, tables[table_name])
		
		#mappers['results_method'] = mapper(ResultsMethod, tables['results_method'], properties={'call_method': relation(CallMethod),\
		#																				'phenotype_method': relation(PhenotypeMethod),\
		#																				'results_method_type': relation(ResultsMethodType)})
	
	_setup_mappers = classmethod(_setup_mappers)


if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = GenomeDatabase
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
