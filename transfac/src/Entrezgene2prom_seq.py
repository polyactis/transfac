#!/usr/bin/env python
"""
Usage: Entrezgene2prom_seq.py [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, sequence(default)
	-t ...,	prom_seq table, 'transfac.prom_seq'(default)
	-b,	enable debug
	-r,	enable report
	-c,	commit this database transaction
	-h, --help	show this help

Examples:
	Entrezgene2prom_seq.py -t transfac.test_prom_seq -c -b -r

Description:
	This program reads the location data of genes in entrezgene_mapping_table,
	get potential regulatory sequence(upstream 10kb and 1st intron) from
	annot_assembly_table(raw_sequence_table) and submit them to prom_seq_table.
	
	NOTE: upstream 10kb will be truncated if there's any other genes closeby. Another
	way to handle this is to mask these coding region with 'N'.
	
	
"""

import psycopg, sys, getopt, os, csv, cStringIO, re
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
sys.path += [os.path.join(os.path.expanduser('~/script/microarray/bin'))]
from codense.common import db_connect, org_short2long, tax_id2org, \
	get_sequence_segment, pg_1d_array2python_ls

class prom_seq_attr:
	def __init__(self):
		self.prom_acc = ''
		self.chromosome = ''
		self.strand = ''
		self.prom_genome_start = 0
		self.prom_genome_end = 0
		self.organism = ''
		self.prom_type_id = 1
		self.sequence = ''
		self.comment = ''

class Entrezgene2prom_seq:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema='sequence',\
		prom_seq_table='transfac.prom_seq', commit=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.prom_seq_table = prom_seq_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def submit_to_prom_seq(self, curs, prom_seq_table, ps_attr_instance):
		curs.execute("insert into %s(prom_acc, chromosome, strand, prom_genome_start, prom_genome_end, \
			organism, prom_type_id, sequence) VALUES('%s', '%s', '%s', %s, %s, \
			'%s', %s, '%s')"%(prom_seq_table, \
			ps_attr_instance.prom_acc, ps_attr_instance.chromosome, ps_attr_instance.strand, ps_attr_instance.prom_genome_start, ps_attr_instance.prom_genome_end,\
			ps_attr_instance.organism, ps_attr_instance.prom_type_id, ps_attr_instance.sequence))
	
	def return_closest_anchor(self, curs, start_or_stop, putative_region, gene_id, tax_id, genomic_gi, \
		entrezgene_mapping_table='entrezgene_mapping'):
		"""
		11-13-05
			For plus stand, we return the maximum of 'stop'
			For minus strand, we return the minimum of 'start'
		"""
		if start_or_stop == 'stop':
			anchor_ls = [putative_region[0]]
		elif start_or_stop == 'start':
			anchor_ls = [putative_region[1]]
		else:
			sys.stderr.write("start_or_stop: %s not supported.\n"%start_or_stop)
			sys.exit(2)
		#same tax_id, same genomic_gi
		curs.execute("select gene_id, %s from %s where  tax_id=%s and genomic_gi=%s and %s>=%s and %s<=%s"%\
			(start_or_stop, entrezgene_mapping_table, tax_id, genomic_gi, start_or_stop, putative_region[0], \
			start_or_stop, putative_region[1]))
		rows = curs.fetchall()
		for row in rows:
			new_gene_id, anchor_pos = row
			if start_or_stop == 'stop':
				anchor_pos += 1
			elif start_or_stop == 'start':
				anchor_pos -= 1
			if new_gene_id!=gene_id:	#a different gene
				anchor_ls.append(anchor_pos)
				if self.debug:
					print "gene_id",gene_id,"start_or_stop",start_or_stop,"putative_region",putative_region
					print "new_gene_id",new_gene_id,"anchor_pos",anchor_pos
		if start_or_stop == 'stop':
			return max(anchor_ls)
		elif start_or_stop == 'start':
			return min(anchor_ls)
		
	def get_prom_seq_from_entrezgene_mapping_table(self, curs, prom_seq_table, entrezgene_mapping_table='entrezgene_mapping', \
		annot_assembly_table = 'annot_assembly'):
		sys.stderr.write("Getting prom_seq from entrezgene_mapping_table...\n")
		curs.execute("DECLARE crs CURSOR FOR SELECT e.gene_id, e.genomic_gi, e.tax_id, a.chromosome, e.strand,\
			e.start, e.stop, e.mrna_start, e.mrna_stop, e.cds_start, e.cds_stop from %s e, %s a \
			where e.genomic_gi=a.gi"%(entrezgene_mapping_table, annot_assembly_table))
		curs.execute("fetch 10000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				gene_id, genomic_gi, tax_id, chromosome, strand, start, stop, mrna_start, mrna_stop, cds_start, cds_stop = row
				seg_loc_ls = []
				if cds_start and cds_stop:
					cds_start = pg_1d_array2python_ls(cds_start, int)
					cds_stop = pg_1d_array2python_ls(cds_stop, int)
					for i in range(len(cds_start)):
						seg_loc_ls.append([cds_start[i],cds_stop[i]])
				elif mrna_start and mrna_stop:
					mrna_start = pg_1d_array2python_ls(mrna_start, int)
					mrna_stop = pg_1d_array2python_ls(mrna_stop, int)
					for i in range(len(mrna_start)):
						seg_loc_ls.append([mrna_start[i],mrna_stop[i]])
				else:
					seg_loc_ls.append([start, stop])
				seg_loc_ls.sort()	#some genes have reversed cds order
				ps_attr_instance = prom_seq_attr()
				ps_attr_instance.prom_acc = gene_id
				ps_attr_instance.chromosome = chromosome
				ps_attr_instance.organism = tax_id2org(tax_id)
				upstream_loc_ls = [0,0]
				instron_1st_loc_ls = []
				if strand=='1':	#plus strand
					ps_attr_instance.strand = '+'
					upstream_loc_ls[1] = seg_loc_ls[0][0]-1
					upstream_loc_ls[0] = upstream_loc_ls[1] - 9999
					if upstream_loc_ls[0]<1:	#in case exceed the chromosome boundary
						upstream_loc_ls[0] = 1
					#check whether there's gene upstream
					upstream_loc_ls[0] = self.return_closest_anchor(curs, 'stop', upstream_loc_ls, gene_id, tax_id, genomic_gi, \
						entrezgene_mapping_table)
					if upstream_loc_ls[0]>upstream_loc_ls[1]:	#No upstream
						if self.debug:
							sys.stderr.write("\tgene_id: %s no upstream\n"%gene_id)
						upstream_loc_ls = []
					if len(seg_loc_ls)>1:	#the first intron
						instron_1st_loc_ls.append(seg_loc_ls[0][1]+1)
						instron_1st_loc_ls.append(seg_loc_ls[1][0]-1)
				elif strand=='-1':	#minus strand
					ps_attr_instance.strand = '-'
					upstream_loc_ls[0] = seg_loc_ls[-1][1]+1
					upstream_loc_ls[1] = upstream_loc_ls[0] + 9999
					#NOTE: exceeding the chromosome boundary is taken care of by get_sequence_segment()
					#check whether there's gene upstream
					upstream_loc_ls[1] = self.return_closest_anchor(curs, 'start', upstream_loc_ls, gene_id, tax_id, genomic_gi, \
						entrezgene_mapping_table)
					if upstream_loc_ls[0]>upstream_loc_ls[1]:	#No upstream
						sys.stderr.write("\tgene_id: %s no upstream\n"%gene_id)
						upstream_loc_ls = []
					if len(seg_loc_ls)>1:	#the first intron
						instron_1st_loc_ls.append(seg_loc_ls[-2][1]+1)
						instron_1st_loc_ls.append(seg_loc_ls[-1][0]-1)
				else:	#ignore genes with no strand info, some are not real genes
					continue
				
				#1st deal with upstream_loc_ls
				if upstream_loc_ls:
					ps_attr_instance.prom_genome_start = upstream_loc_ls[0]
					ps_attr_instance.prom_genome_end = upstream_loc_ls[1]
					ps_attr_instance.prom_type_id = 1
					ps_attr_instance.sequence = get_sequence_segment(curs, genomic_gi, upstream_loc_ls[0], upstream_loc_ls[1])
					self.submit_to_prom_seq(curs, prom_seq_table, ps_attr_instance)
				#2nd handle instron_1st_loc_ls, might not exist
				if instron_1st_loc_ls:
					if instron_1st_loc_ls[0]>instron_1st_loc_ls[1]:
						sys.stderr.write("\tgene_id: %s weird 1st intron %s.\n"%(gene_id, instron_1st_loc_ls))
					ps_attr_instance.prom_genome_start = instron_1st_loc_ls[0]
					ps_attr_instance.prom_genome_end = instron_1st_loc_ls[1]
					ps_attr_instance.prom_type_id = 5
					ps_attr_instance.sequence = get_sequence_segment(curs, genomic_gi, instron_1st_loc_ls[0], instron_1st_loc_ls[1])
					self.submit_to_prom_seq(curs, prom_seq_table, ps_attr_instance)
				counter += 1
			if self.report:
				sys.stderr.write("%s\t%s"%('\x08'*20, counter))
			if self.debug:	#enough
				break
			curs.execute("fetch 10000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done getting prom_seq from entrezgene_mapping_table.\n")
		
	def run(self):
		"""
		11-14-05
			--db_connect()
			--get_prom_seq_from_entrezgene_mapping_table()
				--return_closest_anchor()
				--tax_id2org()
				--get_sequence_segment()
				--submit_to_prom_seq()
		"""
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		self.get_prom_seq_from_entrezgene_mapping_table(curs, self.prom_seq_table)
		if self.commit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:cbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'sequence'
	prom_seq_table = 'transfac.prom_seq'
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
			prom_seq_table = arg
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	instance = Entrezgene2prom_seq(hostname, dbname, schema, prom_seq_table, \
		commit, debug, report)
	instance.run()
