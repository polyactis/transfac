#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	GenomeDB.py -u crocea -k genome
	
	#setup database in mysql
	GenomeDB.py -v mysql -u yh -z papaya -d genome -k ""
	
Description:
	2008-07-09
	This is a wrapper for the genome database, build on top of elixir. supercedes the table definitions in genomedb.sql.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule.db.GenomeDB import *

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = GenomeDatabase
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.setup()
	
	#2008-10-01	get gene model and pickle it into a file
	gene_id2model, chr_id2gene_id_ls = instance.get_gene_id2model()
	from pymodule import PassingData
	gene_annotation = PassingData()
	gene_annotation.gene_id2model = gene_id2model
	gene_annotation.chr_id2gene_id_ls = chr_id2gene_id_ls
	import cPickle
	picklef = open('/tmp/at_gene_model_pickelf', 'w')
	cPickle.dump(gene_annotation, picklef, -1)
	picklef.close()
	
	import sqlalchemy as sql
	#print dir(Gene)
	#print Gene.table.c.keys()
	#results = instance.session.query(Gene)
	#results = instance.session.execute(sql.select([Gene.table]))
	#print dir(results)
	#print dir(Gene.query)
	
	import pdb
	pdb.set_trace()
	i = 0
	block_size = 10
	rows = Gene.query.offset(i).limit(block_size)
	print dir(rows)
	while rows.count()!=0:
		print rows.count()
		for row in rows:
			i += 1
			print row.gene_id, row.gene_symbol
		if i>=5*block_size:
			break
		rows = Gene.query.offset(i).limit(block_size)
