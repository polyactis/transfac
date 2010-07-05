#!/usr/bin/env python
"""
Usage: AnalyzeTRANSFACHits.py -i input_dir -o output_prefix [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, transfac(default)
	-i ...,	input_dir
	-o ...,	output_prefix
	-m ..., 	matrix2no_of_random_hits_table, 'transfac.matrix2no_of_random_hits'(default)
	-c	commit the database transaction
	-b	enable debugging, no debug by default
	-r	report the progress(a number)
	-h, --help              show this help

Examples:
	AnalyzeTRANSFACHits.py -i tmp/tmp -o tmp/test_AnalyzeTRANSFACHits

Description:
	A program to calculate empirical p-values for TRANSFAC hits.
	Further analyze the FDR(Storey2003).
	~/pickle/mt_id_gc_perc2no_of_random_hits.pickle is a random hits files
"""


import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv, cPickle, rpy
from rpy import r
from codense.common import db_connect

import matplotlib as mpl; mpl.use("Agg")
from matplotlib import rcParams
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 4
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 4
rcParams['axes.titlesize'] = 6
rcParams['xtick.labelsize'] = 4
rcParams['ytick.labelsize'] = 4
rcParams['lines.linewidth'] = 0.5	#2008-10-31 make the linewidth of boxplot smaller

import pylab


class AnalyzeTRANSFACHits:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='transfac', 
		input_dir=None, output_prefix=None, matrix2no_of_random_hits_table='transfac.matrix2no_of_random_hits',\
		commit=0, debug=0, report=0):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_dir = input_dir
		self.output_prefix = output_prefix
		self.matrix2no_of_random_hits_table = matrix2no_of_random_hits_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
		
		self.GC_dict = {'G':1, 'C':1, 'g':1, 'c':1}
		self.p_value_list = []	#to store all p-values
		
	def get_mt_id_gc_perc2no_of_random_hits(self, curs, \
		matrix2no_of_random_hits_table='transfac.matrix2no_of_random_hits'):
		
		sys.stderr.write("Getting mt_id_gc_perc2no_of_random_hits...\n")
		curs.execute("DECLARE crs CURSOR FOR select mt_id, gc_perc, no_of_hits\
			from %s"%matrix2no_of_random_hits_table)
		curs.execute("fetch 10000 from crs")
		rows = curs.fetchall()
		mt_id_gc_perc2no_of_random_hits = {}
		counter = 0
		while rows:
			for row in rows:
				mt_id, gc_perc, no_of_hits = row
				gc_perc = int(gc_perc*10)	#round to integer, float is hard to do dictionary match
				key = (mt_id, gc_perc)
				if key not in mt_id_gc_perc2no_of_random_hits:
					mt_id_gc_perc2no_of_random_hits[key] = []
				mt_id_gc_perc2no_of_random_hits[key].append(no_of_hits)
				counter += 1
			if counter%5000 == 0:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
			curs.execute("fetch 10000 from crs")
			rows = curs.fetchall()
		#sort it
		for mt_id in mt_id_gc_perc2no_of_random_hits:
			mt_id_gc_perc2no_of_random_hits[mt_id].sort()
		
		sys.stderr.write("Done.\n")
		return mt_id_gc_perc2no_of_random_hits
	
	def get_seq_id_gc_percentage_length(self, curs, seq_id, prom_seq_table='transfac.prom_seq'):
		"""
		01-18-06
			global structure used:
				self.GC_dict
		"""
		curs.execute("select sequence from %s where id=%s"%(prom_seq_table, seq_id))
		rows = curs.fetchall()
		sequence = rows[0][0]
		no_of_GCs = 0.0
		for letter in sequence:
			if letter in self.GC_dict:
				no_of_GCs += 1
		if len(sequence):	#could be 0
			return no_of_GCs/len(sequence), len(sequence)
		else:
			return None, None
	
	def get_hit_pvalue(self, mt_id, gc_perc, no_of_hits, sequence_length, mt_id_gc_perc2no_of_random_hits):
		"""
		01-18-06
			1. mt_id_gc_perc2no_of_random_hits is based on 10kb sequence
				so need adjustment
			2. if (mt_id, gc_perc) doesn't appear in mt_id_gc_perc2no_of_random_hits, return None
		
		global structure used:
			self.log_f
		"""
		pvalue = None
		if no_of_hits == 0:	#mostly it won't happen, 
			pvalue = 1
		else:
			adjusted_no_of_hits = no_of_hits*10000.0/sequence_length
			gc_perc_to_refer = int(gc_perc*10)	#round to integer
			key = (mt_id, gc_perc_to_refer)
			if key not in mt_id_gc_perc2no_of_random_hits:
				self.log_f.write("\t Warning: %s doesn't appear in mt_id_gc_perc2no_of_random_hits.\n"%repr(key))
			else:
				no_of_random_hits_list = mt_id_gc_perc2no_of_random_hits[key]
				length_no_of_random_hits_list = len(no_of_random_hits_list)
				for i in range(length_no_of_random_hits_list):
					if no_of_random_hits_list[i]>=adjusted_no_of_hits:
						break
				if i==length_no_of_random_hits_list-1 and no_of_random_hits_list[i]<adjusted_no_of_hits:
					pvalue = 0	#adjusted_no_of_hits is larger than any in no_of_random_hits_list
				else:
					pvalue = (length_no_of_random_hits_list-i)/10000.0	#10k random sequences total
						#length_no_of_random_hits_list might <10000, because no_of_random_hits_list 
						#doesn't include 0 hits.
		return pvalue

	def write_down_mt_id2no_of_hits(self, curs, writer, seq_id, mt_id2no_of_hits,\
		mt_id_gc_perc2no_of_random_hits):
		
		gc_perc, sequence_length = self.get_seq_id_gc_percentage_length(curs, seq_id)
		if sequence_length:	#not None
			for mt_id, no_of_hits in mt_id2no_of_hits.iteritems():
				pvalue = self.get_hit_pvalue(mt_id, gc_perc, no_of_hits, sequence_length, mt_id_gc_perc2no_of_random_hits)
				if pvalue:	#could be None
					writer.writerow([seq_id, sequence_length, mt_id, gc_perc, no_of_hits, pvalue])
					self.p_value_list.append(pvalue)
		
	def parse_file(self, curs, input_fname, writer, mt_id_gc_perc2no_of_random_hits):
		sys.stderr.write("\t Parsing %s..."%os.path.basename(input_fname))
		inf = open(input_fname, 'r')
		mt_id2no_of_hits = {}
		seq_id = None
		for line in inf:
			if line[:10]=='Inspecting':
				if mt_id2no_of_hits and seq_id:	#1st time, both mt_id2no_of_hits and seq_id are nothing
					self.write_down_mt_id2no_of_hits(curs, writer, seq_id, mt_id2no_of_hits, mt_id_gc_perc2no_of_random_hits)
				seq_id = int(line.split()[-1])	#\n is discarded automatically by split()
				mt_id2no_of_hits = {}	#clear the dictionary
			elif line[0] == ' ' and line[:6]!=' Total' and line[:6]!=' Frequ':	#the first character is blank, but exclude the end statistic part	
				ls = line[:-1].split('|')
				mt_id = ls[0].strip()	#remove spaces
				"""
				bs_disp_start_strand = ls[1].strip()
				bs_disp_start = int(bs_disp_start_strand[:-3])
				strand = bs_disp_start_strand[-2]
				core_similarity_score = ls[2]	#01-07-06
				matrix_similarity_score = ls[3]	#01-07-06
				sequence = ls[4].strip()
				#01-03-06
				bs_disp_end = bs_disp_start + len(sequence) - 1
				"""
				if mt_id not in mt_id2no_of_hits:
					mt_id2no_of_hits[mt_id] = 0
				mt_id2no_of_hits[mt_id] += 1
		
		if mt_id2no_of_hits and seq_id:	#don't forget the last sequence
			self.write_down_mt_id2no_of_hits(curs, writer, seq_id, mt_id2no_of_hits, mt_id_gc_perc2no_of_random_hits)
		del inf
		sys.stderr.write("Done.\n")
	
	def remove_top_p_values(self, p_value_list, top_p_value_cutoff):
		"""
		01-18-06
			p_value_list already sorted
			throw away top bad p_value's and fill in proportional amount
		"""
		while p_value_list[-1]>top_p_value_cutoff:	#starting from the end
			p_value_list.pop()
		#fill in same amount of fake values to lower the p_value histogram tail
		lower_p_value_cutoff = top_p_value_cutoff-(1-top_p_value_cutoff)
		for i in range(len(p_value_list)-1, -1, -1):	#reverse
			if p_value_list[i]>lower_p_value_cutoff:
				p_value_list.append(top_p_value_cutoff+(1-top_p_value_cutoff)/2)
			else:
				break
		
	def draw_pvalue_histogram(self, p_value_list, figure_fname):
		"""
		2008-11-30
			use pylab to draw histograms
		"""
		sys.stderr.write("Drawing p_value_list histogram...")
		#r.png(figure_fname)
		#r.hist(p_value_list, main='histogram',xlab='p_value',ylab='freq')
		#r.dev_off()
		pylab.clf()
		pylab.hist(p_value_list, 100, alpha=0.4)
		pylab.title('histogram of pvalue')
		pylab.xlabel('p_value')
		pylab.savefig(figure_fname, dpi=300)
		sys.stderr.write("Done.\n")
	
	def calculate_pi0_list(self, p_value_list, figure_fname, top_p_value_cutoff=1, lambda_gap=0.001):
		"""
		2008-11-30
			add option lambda_gap
			use pylab to draw plots
		01-18-06
			Details see Remark B in Storey2003 PNAS.
			
			fomula:
				pi0 = #{p>lambda}/(m(1-lambda))
			
			top_p_value_cutoff is used to cut off many too high p-values in remove_top_p_values()
		"""
		sys.stderr.write("Calculating pi0 (proportion of truly null features) ...\n")
		#find all lambda's and corresponding #{p>lambda}
		lambda_to_cmp = 0
		lambda_gap = lambda_gap
		lambda_list = []
		no_of_pvalues_above_lambda_list = []
		no_of_total_p_values = len(p_value_list)
		for i in range(no_of_total_p_values):
			if p_value_list[i]>=lambda_to_cmp:
				lambda_list.append(lambda_to_cmp)
				no_of_pvalues_above_lambda_list.append(no_of_total_p_values-i)
				lambda_to_cmp += lambda_gap
				if lambda_to_cmp>=top_p_value_cutoff:
					break
		#calculate pi0
		pi0_list = []
		for i in range(len(lambda_list)):
			pi0 = no_of_pvalues_above_lambda_list[i]/(no_of_total_p_values*(1-lambda_list[i]))	#lambda_list is float
			pi0_list.append(pi0)
		
		print "\t lambda_list", lambda_list
		print "\t no_of_pvalues_above_lambda_list", no_of_pvalues_above_lambda_list
		print "\t pi0_list", pi0_list
		#r.png(figure_fname)
		#r.plot(lambda_list, pi0_list, xlab='lambda', ylab='pi0')
		#r.dev_off()
		pylab.clf()
		pylab.plot(lambda_list, pi0_list, 'o-', **self.plot_kw)
		pylab.title(r'estimate $\pi_0$')
		pylab.xlabel(r'$\lambda$')
		pylab.ylabel(r'$\pi_0$')
		pylab.savefig(figure_fname, dpi=300)
		sys.stderr.write("Done.\n")
		
		return lambda_list, pi0_list
	
	def estimate_pi0(self, lambda_list, pi0_list):
		"""
		01-19-06
			Storey2003, (natural) cubic spline, df=3
		"""
		sys.stderr.write("Estimating pi0...\n")
		rpy.set_default_mode(rpy.NO_CONVERSION)
		s = r.smooth_spline(lambda_list, pi0_list, df=3)
		rpy.set_default_mode(rpy.BASIC_CONVERSION)
		estimated_pi0 = r.predict(s,1)['y']
		print "\t estimated_pi0:", estimated_pi0
		sys.stderr.write("Done.\n")
		return estimated_pi0
	
	def cal_q_value_list(self, p_value_list, estimated_pi0, top_p_value_cutoff, output_prefix):
		"""
		2008-11-30
			use pylab to draw plots
		01-19-06
		"""
		sys.stderr.write("Calculating q_value list...\n")
		used_p_value_list = []
		q_value_list = []
		no_of_sig_hits_list = []
		no_of_total_p_values = len(p_value_list)
		prev_q_value = None
		for i in range(no_of_total_p_values-1, -1, -1):
			if p_value_list[i] < top_p_value_cutoff:
				current_q_value = estimated_pi0*no_of_total_p_values*p_value_list[i]/(i+1)
				if prev_q_value:
					current_q_value = min(current_q_value, prev_q_value)
					prev_q_value = current_q_value
				else:
					prev_q_value = current_q_value
				q_value_list.append(current_q_value)
				used_p_value_list.append(p_value_list[i])
				no_of_sig_hits_list.append(i+1)
		#print "\t q_value_list:", q_value_list
		#print "\t used_p_value_list:", used_p_value_list
		#print "\t no_of_sig_hits_list:", no_of_sig_hits_list
		
		figure_fname = '%s.p_value.vs.q_value.png'%output_prefix
		#r.png(figure_fname)
		#r.plot(q_value_list, used_p_value_list, main='p_value.vs.q_value', xlab='q_value', ylab='p_value')
		#r.dev_off()
		pylab.clf()
		pylab.plot(q_value_list, used_p_value_list, 'o', **self.plot_kw)
		pylab.title('p_value vs q_value')
		pylab.xlabel('q_value')
		pylab.ylabel('p_value')
		pylab.savefig(figure_fname, dpi=300)
		
		for n in [10000, 5000, 1000, 500, 100]:
			figure_fname = '%s.p_value.vs.q_value_zoom_%s.png'%(output_prefix,n)
			pylab.clf()
			pylab.plot(q_value_list[-n:], used_p_value_list[-n:], 'o', **self.plot_kw)
			pylab.title('p_value vs q_value')
			pylab.xlabel('q_value')
			pylab.ylabel('p_value')
			pylab.savefig(figure_fname, dpi=300)
		
		figure_fname = '%s.no_of_sig_hits.vs.q_value.png'%output_prefix
		#r.png(figure_fname)
		#r.plot(q_value_list, no_of_sig_hits_list, main='no_of_sig_hits.vs.q_value', xlab='q_value', ylab='no_of_sig_hits')
		#r.dev_off()
		pylab.clf()
		pylab.plot(q_value_list, no_of_sig_hits_list, 'o', **self.plot_kw)
		pylab.title('no_of_sig_hits vs q_value')
		pylab.xlabel('q_value')
		pylab.ylabel('no_of_sig_hits')
		pylab.savefig(figure_fname, dpi=300)
		
		#2008-11-30	draw close-up of the plot above
		for n in [10000, 5000, 1000, 500, 100]:
			figure_fname = '%s.no_of_sig_hits.vs.q_value_zoom_%s.png'%(output_prefix, n)
			pylab.clf()
			pylab.plot(q_value_list[-n:], no_of_sig_hits_list[-n:], 'o', **self.plot_kw)
			pylab.title('no_of_sig_hits vs q_value')
			pylab.xlabel('q_value')
			pylab.ylabel('no_of_sig_hits')
			pylab.savefig(figure_fname, dpi=300)
		
		sys.stderr.write("Done.\n")
	
	plot_kw = {'linewidth':0.5,\
				'markerfacecolor': 'w',\
				'markersize': 1,\
				'alpha':0.6}
	
	def run(self):
		"""
		01-18-06
		
			--db_connect()
			--get_mt_id_gc_perc2no_of_random_hits()
			--parse_file()
				--write_down_mt_id2no_of_hits()
					--get_seq_id_gc_percentage_length()
					--get_hit_pvalue()
			
			--draw_pvalue_histogram()
			
			--calculate_pi0()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		data_fname = '%s.data'%self.output_prefix
		if os.path.isfile(data_fname):
			sys.stderr.write("Getting p_value from %s..."%os.path.basename(data_fname))
			reader = csv.reader(open(data_fname), delimiter='\t')
			for row in reader:
				self.p_value_list.append(float(row[5]))
			del reader
			sys.stderr.write("Done.\n")
		else:
			pickle_fname = os.path.expanduser('~/pickle/mt_id_gc_perc2no_of_random_hits.pickle')
			if os.path.isfile(pickle_fname):
				mt_id_gc_perc2no_of_random_hits = cPickle.load(open(pickle_fname))
			else:
				mt_id_gc_perc2no_of_random_hits = self.get_mt_id_gc_perc2no_of_random_hits(curs,\
					self.matrix2no_of_random_hits_table)
				of = open(pickle_fname, 'w')
				cPickle.dump(mt_id_gc_perc2no_of_random_hits, of)
				del of
			writer = csv.writer(open(data_fname, 'w') , delimiter='\t')
			
			self.log_f = open('%s.log'%self.output_prefix,'w')
			
			files = os.listdir(self.input_dir)
			files.sort()
			sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
			for input_fname in files:
				input_fname = os.path.join(self.input_dir, input_fname)
				self.parse_file(curs, input_fname, writer, mt_id_gc_perc2no_of_random_hits)
			del writer
			self.log_f.close()
		
		self.p_value_list.sort()
		top_p_value_cutoff = 0.95	#important not 1, p_value histogram shows an abnormal peak from 0.95 to 1
		top_p_value_list = self.remove_top_p_values(self.p_value_list, top_p_value_cutoff)
		
		figure_fname = '%s_p_value_hist.png'%self.output_prefix
		self.draw_pvalue_histogram(self.p_value_list, figure_fname)
		
		figure_fname = '%s_pi0Tolambda.png'%self.output_prefix
		lambda_list, pi0_list = self.calculate_pi0_list(self.p_value_list, figure_fname, top_p_value_cutoff)
		
		estimated_pi0 = self.estimate_pi0(lambda_list, pi0_list)
		
		self.cal_q_value_list(self.p_value_list, estimated_pi0, top_p_value_cutoff, self.output_prefix)
		
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:m:cbr", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'transfac'
	input_dir = None
	output_prefix = None
	matrix2no_of_random_hits_table = 'transfac.matrix2no_of_random_hits'
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
		elif opt in ("-i",):
			input_dir = arg
		elif opt in ("-o",):
			output_prefix = arg
		elif opt in ("-m",):
			matrix2no_of_random_hits_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
		
	if input_dir and output_prefix:
		instance = AnalyzeTRANSFACHits(hostname, dbname, schema, input_dir, output_prefix, \
			matrix2no_of_random_hits_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
