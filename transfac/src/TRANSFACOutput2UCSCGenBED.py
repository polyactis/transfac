#!/usr/bin/env python
"""
Usage: TRANSFACOutput2UCSCGenBED.py -i input_dir -o output_fname [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, sequence(default)
	-i ...,	input_dir, where TRANSFAC output files are
	-o ...,	output_fname
	-x ...,	tax_id, 9606(human, default)
	-b	enable debugging, no debug by default
	-r	report the progress(a number)
	-h, --help              show this help

Examples:
	~/script/transfac/src/TRANSFACOutput2UCSCGenBED.py -i transfac/genome_seq_hs_out/
		-o transfac/genome_seq_hs_out.UCSCGen.BED
Description:
	Program to collect all output files of TRANSFAC's MATCH, reformat them
	into BED format of UCSC Genome Browser for custom track.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import os, sys, getopt, random, csv
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
from codense.common import db_connect, get_gene_id2gene_symbol

class TRANSFACOutput2UCSCGenBED:
	"""
	2006-12-14
		
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='sequence', input_dir=None, \
		output_fname=None, tax_id=9606, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_dir = input_dir
		self.output_fname = output_fname
		self.tax_id = int(tax_id)
		self.debug = int(debug)
		self.report = int(report)
		
		#based on HTML color names from http://en.wikipedia.org/wiki/Web_colors, except white
		self.color_code_list = ['0,255,255', '0,0,0', '0,0,255', '255,0,255', \
			'128,128,128', '0,128,0', '0,255,0', '128,0,0',\
			'0,0,128', '128,128,0', '128,0,128', '255,0,0',\
			'192,192,192', '0,128,128', '255,255,0']
	
	def get_id2chr_start_stop(self, curs, tax_id, raw_sequence_table='raw_sequence', annot_assembly_table='annot_assembly'):
		"""
		2006-12-14
			get the mapping between prom_id and chromosome, start, stop
		"""
		sys.stderr.write("Getting id2chr_start_stop...\n")
		curs.execute("DECLARE crs CURSOR FOR SELECT r.id, a.chromosome, r.start, r.stop from %s r, %s a\
				where r.acc_ver=a.acc_ver and a.tax_id=%s"%(raw_sequence_table, annot_assembly_table, tax_id))
		curs.execute("fetch 2000 from crs")
		rows = curs.fetchall()
		counter = 0
		id2chr_start_stop = {}
		while rows:
			for row in rows:
				id, chromosome, start, stop = row
				id2chr_start_stop[id] = ['chr%s'%chromosome, start, stop]
			counter += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			curs.execute("fetch 2000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("done.\n")
		return id2chr_start_stop
	
	def get_mt_id2gene_symbol(self, curs, gene_id2gene_symbol, tax_id, table='transfac.mt2gene_id'):
		"""
		2006-12-14
			To avoid different color assigned to the same transfactor with multi-matrix-ids, connect mt_id to gene_symbol
		"""
		sys.stderr.write("Getting mt_id2gene_symbol...\n")
		curs.execute("select mt_id, gene_id from %s where tax_id=%s"%(table, tax_id))
		rows = curs.fetchall()
		mt_id2gene_symbol = {}
		for row in rows:
			mt_id, gene_id = row
			mt_id2gene_symbol[mt_id] = gene_id2gene_symbol[gene_id]
		sys.stderr.write("done.\n")
		return mt_id2gene_symbol
	
	def get_mt_id2color_code(self, input_dir, color_code_list, mt_id2gene_symbol):
		"""
		2006-12-14
			assign the color_code to different mt_id ordered by frequency they appear
		"""
		sys.stderr.write("Getting tfbs2color_code...\n")
		files = os.listdir(input_dir)
		files.sort()
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		mt_id_gene_symbol2freq = {}
		for i in range(len(files)):
			f = files[i]
			sys.stderr.write('\t %s/%s %s'%((i+1), len(files), f))
			f_path = os.path.join(input_dir, f)
			inf = open(f_path, 'r')
			for line in inf:
				if line[:10]=='Inspecting':
					prom_id = int(line.split()[-1])	#\n is discarded automatically by split()
				elif line[0] == ' ' and line[:6]!=' Total' and line[:6]!=' Frequ':	#the first character is blank, but exclude the end statistic part	
					ls = line[:-1].split('|')
					mt_id = ls[0].strip()	#remove spaces
					if mt_id in mt_id2gene_symbol:
						mt_id_gene_symbol = mt_id2gene_symbol[mt_id]	#get the gene symbol
						if mt_id_gene_symbol not in mt_id_gene_symbol2freq:
							mt_id_gene_symbol2freq[mt_id_gene_symbol] = 0
						mt_id_gene_symbol2freq[mt_id_gene_symbol] += 1
			del inf
			sys.stderr.write(".\n")
		freq_mt_id_gene_symbol_ls = []
		for mt_id_gene_symbol, freq in mt_id_gene_symbol2freq.iteritems():
			freq_mt_id_gene_symbol_ls.append((freq, mt_id_gene_symbol))
		freq_mt_id_gene_symbol_ls.sort()
		
		mt_id_gene_symbol2color_code = {}
		for i in range(len(freq_mt_id_gene_symbol_ls)):
			mt_id_gene_symbol2color_code[freq_mt_id_gene_symbol_ls[i][1]] = color_code_list[i%len(color_code_list)]	#recycle the color_code_list
		
		sys.stderr.write("Done.\n")
		return mt_id_gene_symbol2color_code
	
	def parse_files(self, input_dir, output_fname, id2chr_start_stop, mt_id_gene_symbol2color_code, mt_id2gene_symbol):
		"""
		2006-12-14
			To avoid different color assigned to the same transfactor with multi-matrix-ids, use gene_symbol instead
		"""
		sys.stderr.write("Parsing files ...\n")
		files = os.listdir(input_dir)
		files.sort()
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		outf = open(output_fname, 'w')
		outf.write('browser position chr1:1000-4000\n')
		outf.write('browser squish knownGene\n')	#track name is the name of the primary table on which the the track is based. To identify this table, open up the Table Browser, select the correct genome assembly, then select the track name from the track list. The table list will show the primary table.
		outf.write('browser squish refGene\n')
		outf.write('browser squish xenoRefGene\n')
		outf.write('browser squish ensGene\n')
		outf.write('browser squish wgRna\n')
		outf.write('browser squish evofold\n')
		outf.write('browser squish all_mrna\n')
		outf.write('browser squish uniGene_2\n')
		outf.write('browser squish altGraphX\n')
		outf.write('browser squish phastCons17way\n')
		outf.write('browser dense rmsk\n')	#repeat masker only for dense
		outf.write('track name="TRANSFAC hits" description="TRANSFAC hits, Repeat Masked in advance" visibility=4 itemRgb="On"\n')
		writer = csv.writer(outf, delimiter = '\t')
		for i in range(len(files)):
			f = files[i]
			sys.stderr.write('\t %s/%s %s'%((i+1), len(files), f))
			f_path = os.path.join(input_dir, f)
			inf = open(f_path, 'r')
			for line in inf:
				if line[:10]=='Inspecting':
					prom_id = int(line.split()[-1])	#\n is discarded automatically by split()
					chromosome, chr_start, chr_stop = id2chr_start_stop[prom_id]
				elif line[0] == ' ' and line[:6]!=' Total' and line[:6]!=' Frequ':	#the first character is blank, but exclude the end statistic part	
					ls = line[:-1].split('|')
					mt_id = ls[0].strip()	#remove spaces
					bs_disp_start_strand = ls[1].strip()
					bs_disp_start = int(bs_disp_start_strand[:-3])
					strand = bs_disp_start_strand[-2]
					core_similarity_score = ls[2]
					matrix_similarity_score = ls[3]
					sequence = ls[4].strip()
					bs_disp_end = bs_disp_start + len(sequence) - 1
					bs_chr_start = chr_start + bs_disp_start - 1
					bs_chr_stop = bs_chr_start + len(sequence) - 1
					if mt_id in mt_id2gene_symbol:
						row = [chromosome, bs_chr_start, bs_chr_stop, mt_id2gene_symbol[mt_id], matrix_similarity_score, strand, bs_chr_start, bs_chr_stop, mt_id_gene_symbol2color_code[mt_id2gene_symbol[mt_id]] ]
						writer.writerow(row)
			del inf
			sys.stderr.write(".\n")
		del outf
		del writer
		sys.stderr.write(".\n")
	
	def run(self):
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		id2chr_start_stop = self.get_id2chr_start_stop(curs, self.tax_id)
		gene_id2gene_symbol = get_gene_id2gene_symbol(curs, self.tax_id)
		mt_id2gene_symbol = self.get_mt_id2gene_symbol(curs, gene_id2gene_symbol, self.tax_id)
		mt_id_gene_symbol2color_code = self.get_mt_id2color_code(self.input_dir, self.color_code_list, mt_id2gene_symbol)
		self.parse_files(self.input_dir, self.output_fname, id2chr_start_stop, mt_id_gene_symbol2color_code, mt_id2gene_symbol)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:x:br", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'sequence'
	input_dir = None
	output_fname = None
	tax_id = 9606
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
		elif opt in ("-i",):
			input_dir = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-x",):
			tax_id = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r",):
			report = 1
		
	if input_dir and output_fname:
		instance = TRANSFACOutput2UCSCGenBED(hostname, dbname, schema, input_dir, output_fname, \
			tax_id, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
