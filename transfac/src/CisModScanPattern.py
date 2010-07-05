#!/usr/bin/env python
"""
Usage: CisModScanPattern.py -p PATTER_ID -o OUTPUT_DIR [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-i ...,	pattern_table, 'hs_fim_65.pattern_hs_fim_65_m5x65s4l5'(default)
	-p ...,	pattern id
	-t ...,	prom_seq table, 'transfac.prom_seq'(default)
	-o ...,	output directory(all temporary files are here)
	-m ...,	matrix.dat filename, '~/script/transfac/data/matrix.dat'(default)
	-f ...,	profile_filename, '~/script/transfac/data/match_data/vertebrates_m0_8_c0_85_highQual.match.prf'(default)
	-b ...,	path for TRANSFAC's match, '~/script/transfac/bin/match'(default)
	-c ...,	path for CisModScan, '~/bin/CisModScanU'(default)
	-k ...,	maximum number of TFs in the PWM file, 10(default)
	-l ...,	module length, 300(default)
	-x ...,	Expected module occurrent rate (default 0.0005)
	-u,	enable debug
	-r,	enable report
	-h, --help	show this help

Examples:
	CisModScanPattern.py -p 3 -o /tmp/cismod -r -k 15

Description:
	A program wrapping almost everything from TRANSFAC's match to CisModScan
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/transfac/src')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/transfac/src')))
import getopt, csv, cPickle, cStringIO
from transfacdb import match_block_iterator
from sets import Set
from codense.common import db_connect
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
	from python2_3 import *
	
class CisModScanPattern:
	def __init__(self, *argument_list):
		"""
		2006-09-01
		"""
		self.hostname, self.dbname, self.pattern_table, self.pattern_id, self.prom_seq_table, self.output_dir, \
		self.profile_filename, self.matrix_data_path, self.match_bin_path, self.cismodscan_binary_path, \
		self.no_of_tfs, self.mod_length,\
		self.expt_ratio, self.debug, self.report = argument_list
		
	def get_masked_seq(self, curs, gene_id_list, prom_seq_table, output_fname):
		if self.report:
			sys.stderr.write("Getting masked sequences...")
		outf = open(output_fname, 'w')
		for gene_id in gene_id_list:
			curs.execute("select masked_seq from %s where prom_acc='%s' and prom_type_id=1"%(prom_seq_table, gene_id))
			rows  = curs.fetchall()
			outf.write(">%s\n"%gene_id)
			outf.write("%s\n"%rows[0][0])
		del outf
		if self.report:
			sys.stderr.write("done.\n")
	
	def run_transfac(self, input_fname, output_fname, match_bin_path, matrix_path, profile_filename):
		if self.report:
			sys.stderr.write("Running TRANSFAC...\n")
		job_ls = [match_bin_path, matrix_path, input_fname, output_fname, profile_filename]
		pid = os.spawnvp(os.P_WAIT, match_bin_path, job_ls)
		inf = open(output_fname)
		iter = match_block_iterator(inf)
		mt_id2no_of_seqs = {}
		block_no =0
		for block in iter:
			block = cStringIO.StringIO(block)
			line = block.readline()
			prom_id = int(line.split()[-1])	#\n is discarded automatically by split()
			mt_id_set = Set()
			for line in block:
				ls = line[:-1].split('|')
				mt_id = ls[0].strip()
				mt_id_set.add(mt_id)
			for mt_id in mt_id_set:
				if mt_id not in mt_id2no_of_seqs:
					mt_id2no_of_seqs[mt_id] = 0
				mt_id2no_of_seqs[mt_id] += 1
			block_no += 1
			if self.report and block_no%1000==0:
				sys.stderr.write('%s%s'%('\x08'*20, block_no))
		del iter, inf
		if self.report:
			sys.stderr.write("Done\n")
		return mt_id2no_of_seqs
	
	def output_transfac_pwm_cismodscan_format(self, curs, mt_id_list, matrix_table, pwm_fname, pwm_id_mapping_fname):
		"""
		2006-09-01
		"""
		if self.report:
			sys.stderr.write("Output TRANSFAC matrices for CisModScan ... ")
		pwm_f = open(pwm_fname, 'w')
		pwm_id_mapping_f = open(pwm_id_mapping_fname, 'w')
		divid_f= lambda x: x/sum(pwm_row)
		for i in range(len(mt_id_list)):
			mt_id = mt_id_list[i]
			pwm_id_mapping_f.write('%s\t%s\n'%(i+1, mt_id))
			curs.execute("select pwm from %s where mt_id='%s'"%(matrix_table, mt_id))
			rows = curs.fetchall()
			pwm_matrix = rows[0][0][2:-2].split('},{')
			pwm_f.write(">%s\n"%(len(pwm_matrix)))
			for j in range(len(pwm_matrix)):
				pwm_row = pwm_matrix[j].split(',')
				pwm_row = map(float, pwm_row)
				pwm_row = map(divid_f, pwm_row)
				pwm_f.write('%.3f\t%.3f\t%.3f\t%.3f\n'%(pwm_row[0], pwm_row[1], pwm_row[2], pwm_row[3]))
		del pwm_f, pwm_id_mapping_f
		if self.report:
			sys.stderr.write("Done\n")
	
	def run_cismodscan(self, cismodscan_binary_path, seq_fname, pwm_fname, output_fname, no_of_tfs, mod_length, expt_ratio=0.0005):
		"""
		2006-09-04
		"""
		if self.report:
			sys.stderr.write("Running CisModScan...")
		job_ls = [cismodscan_binary_path, '-i', seq_fname, '-o', output_fname, '-M', pwm_fname, '-K', repr(no_of_tfs), '-L', repr(mod_length), '-x', repr(expt_ratio)]
		pid = os.spawnvp(os.P_WAIT, cismodscan_binary_path, job_ls)
		if self.report:
			sys.stderr.write("Done\n")
	
	def get_gene_id_list(self, curs, pattern_table, pattern_id):
		if self.report:
			sys.stderr.write("Getting gene_id_list...")
		curs.execute("select vertex_set from %s where id=%s"%(pattern_table, pattern_id))
		rows = curs.fetchall()
		gene_id_list = rows[0][0][1:-1].split(',')
		gene_id_list = map(int, gene_id_list)
		if self.report:
			sys.stderr.write("Done\n")
		return gene_id_list
	
	def get_top_mt_id_list(self, mt_id2no_of_seqs, no_of_tfs):
		"""
		2006-09-01
			based on how many sequences have this mt_id
		"""
		if self.report:
			sys.stderr.write("Getting top mt_id_list...")
		no_of_seqs_mt_id_list = []
		for mt_id, no_of_seqs in mt_id2no_of_seqs.iteritems():
			no_of_seqs_mt_id_list.append((no_of_seqs, mt_id))
		no_of_seqs_mt_id_list.sort()
		no_of_seqs_mt_id_list.reverse()
		mt_id_list = []
		for i in range(min(no_of_tfs, len(no_of_seqs_mt_id_list))):
			mt_id_list.append(no_of_seqs_mt_id_list[i][1])
		if self.report:
			sys.stderr.write("Done\n")
		return mt_id_list
	
	def run(self):
		"""
		2006-09-04
			-db_connect()
			-get_gene_id_list()
			-get_masked_seq()
			-run_transfac()
			-get_top_mt_id_list()
			-output_transfac_pwm_cismodscan_format()
			-run_cismodscan()
		"""
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		seq_fname = os.path.join(self.output_dir, 'pattern_%s.seq'%(self.pattern_id))
		transfac_output_fname = os.path.join(self.output_dir, 'pattern_%s.match'%(self.pattern_id))
		pwm_fname = os.path.join(self.output_dir, 'pattern_%s.pwm'%(self.pattern_id))
		pwm_id_mapping_fname = os.path.join(self.output_dir, 'pattern_%s.pwm_id_mapping'%(self.pattern_id))
		cismodscan_output_fname = os.path.join(self.output_dir, 'pattern_%s.cismodscan'%(self.pattern_id))
		
		(conn, curs) =  db_connect(self.hostname, self.dbname)
		gene_id_list = self.get_gene_id_list(curs, self.pattern_table, self.pattern_id)
		self.get_masked_seq(curs, gene_id_list, self.prom_seq_table, seq_fname)
		mt_id2no_of_seqs = self.run_transfac(seq_fname, transfac_output_fname, self.match_bin_path, self.matrix_data_path, self.profile_filename)
		mt_id_list = self.get_top_mt_id_list(mt_id2no_of_seqs, self.no_of_tfs)
		matrix_table = 'transfac.matrix'
		self.output_transfac_pwm_cismodscan_format(curs, mt_id_list, matrix_table, pwm_fname, pwm_id_mapping_fname)
		self.run_cismodscan(self.cismodscan_binary_path, seq_fname, pwm_fname, cismodscan_output_fname, self.no_of_tfs, self.mod_length, self.expt_ratio)
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:i:p:t:o:f:m:b:c:k:l:x:ur", ["help", "debug", "report"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	pattern_table = 'hs_fim_65.pattern_hs_fim_65_m5x65s4l5'
	pattern_id = None
	prom_seq_table = 'transfac.prom_seq'
	output_dir = None
	profile_filename = os.path.expanduser('~/script/transfac/data/match_data/vertebrates_m0_8_c0_85_highQual.match.prf')
	matrix_data_path = os.path.expanduser('~/script/transfac/data/matrix.dat')
	match_bin_path = os.path.expanduser('~/script/transfac/bin/match')
	cismodscan_binary_path = os.path.expanduser('~/bin/CisModScanU')
	no_of_tfs = 10
	mod_length = 300
	expt_ratio = 0.0005
	debug = 0
	report = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z",):
			hostname = arg
		elif opt in ("-d",):
			dbname = arg
		elif opt in ("-i", ):
			pattern_table = arg
		elif opt in ("-p", ):
			pattern_id = int(arg)
		elif opt in ("-t", ):
			prom_seq_table = arg
		elif opt in ("-f", ):
			profile_filename = arg
		elif opt in ("-o", ):
			output_dir = arg
		elif opt in ("-m", ):
			matrix_data_path = arg
		elif opt in ("-b", ):
			match_bin_path = arg
		elif opt in ("-c", ):
			cismodscan_binary_path = arg
		elif opt in ("-k", ):
			no_of_tfs = int(arg)
		elif opt in ("-l", ):
			mod_length = int(arg)
		elif opt in ("-x", ):
			expt_ratio = float(arg)
		elif opt in ("-u", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	if pattern_id!=None and output_dir:
		instance = CisModScanPattern(hostname, dbname, pattern_table, pattern_id, prom_seq_table, output_dir, \
			profile_filename, matrix_data_path, match_bin_path, cismodscan_binary_path, no_of_tfs, mod_length, \
			expt_ratio, debug,\
			report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
