#!/usr/bin/env python
"""
Usage: ReformatTRANSFACOutput.py -i input_dir -o output_fname [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database(IGNORE)
	-i ...,	input_dir, where TRANSFAC output files are
	-o ...,	output_fname
	-b	enable debugging, no debug by default
	-r	report the progress(a number)
	-h, --help              show this help

Examples:
	ReformatTRANSFACOutput.py -i ./prom_seq_hs_match_out_m0_8_c0_85/ 
		-o prom_seq_hs_match_out_m0_8_c0_85.reformat

Description:
	Program to collect all output files of TRANSFAC's MATCH, reformat them
	into csv format and append order(based on matrix_similarity_score, 
	core_similarity_score) to each hit.
	
	Two intermediary files: output_fname.int and output_fname.csv
	
"""
import os, sys, getopt, random, csv
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
from codense.common import db_connect


class ReformatTRANSFACOutput:
	"""
	01-06-06
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, input_dir=None, \
		output_fname=None, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_dir = input_dir
		self.output_fname = output_fname
		self.debug = int(debug)
		self.report = int(report)
	
	def Transform2csvFormat(self, input_dir, output_fname):
		"""
		01-07-06
			don't cast matrix_similarity_score and core_similarity_score as float, just treat them as
				0.000 format string(default MATCH output format)
			re-arrange the order of the outputted items in each row
				so that 'sort' could be used to sort by mt_id, matrix_similarity_score, core_similarity_score
			append '-r -n' to 'sort'
		"""
		sys.stderr.write("Transforming TRANSFAC output into csv format...\n")
		files = os.listdir(input_dir)
		files.sort()
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		intermediary_fname = '%s.int'%output_fname
		writer = csv.writer(open(intermediary_fname, 'w'), delimiter = '\t')
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
					bs_disp_start_strand = ls[1].strip()
					bs_disp_start = int(bs_disp_start_strand[:-3])
					strand = bs_disp_start_strand[-2]
					core_similarity_score = ls[2]	#01-07-06
					matrix_similarity_score = ls[3]	#01-07-06
					sequence = ls[4].strip()
					#01-03-06
					bs_disp_end = bs_disp_start + len(sequence) - 1
					row = [mt_id, matrix_similarity_score, core_similarity_score, prom_id, strand, bs_disp_start, bs_disp_end, sequence] #01-07-06
					writer.writerow(row)
			del inf
		del writer
		sys.stderr.write(".\n")
		
		sys.stderr.write("Starting to sort %s..."%intermediary_fname)
		commandline = 'sort -r -n %s > %s'%(intermediary_fname, output_fname)	#01-07-06 reverse and numerical sort
		exit_code = os.system(commandline)
		os.remove(intermediary_fname)
		sys.stderr.write("Done.\n")


	def output_sorted_hits(self, writer, hits_of_one_mt_id):
		"""
		01-07-06
			defunct
		"""
		hits_of_one_mt_id.sort()
		hits_of_one_mt_id.reverse()	#order by matrix_similarity_score, then by core_similarity_score
		for i in range(len(hits_of_one_mt_id)):
			matrix_similarity_score, core_similarity_score, mt_id, prom_id, strand, bs_disp_start, bs_disp_end, sequence = hits_of_one_mt_id[i]
			writer.writerow([mt_id, prom_id, strand, bs_disp_start, bs_disp_end, core_similarity_score, matrix_similarity_score, sequence, i+1])

	def orderReformattedTRANSFACHits(self, ReformatTRANSFACOutput_fname, output_fname):
		"""
		01-06-06
			order the output_fname above based on matrix_similarity_score, then by core_similarity_score
			append the rank to the end of the new row
		01-07-06
			due to memory bottleneck, ordering within one mt_id block is transferred to 'sort' (Transform2csvFormat())
			here just append the rank
		"""
		sys.stderr.write("Append the order of TRANSFAC hits...\n")
		reader = csv.reader(open(ReformatTRANSFACOutput_fname, 'r'), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter ='\t')
		prev_mt_id = None
		counter = 0
		mt_id_counter = 0
		counter_within_one_mt_id_block = 0
		counter_mt_id_ls = []
		for row in reader:
			mt_id, matrix_similarity_score, core_similarity_score, prom_id, strand, bs_disp_start, bs_disp_end, sequence = row
			if prev_mt_id == None:
				prev_mt_id = mt_id
			if mt_id != prev_mt_id:
				counter_mt_id_ls.append([counter_within_one_mt_id_block, prev_mt_id])
				counter_within_one_mt_id_block = 0	#01-07-06 reset the counter_within_one_mt_id_block after above line
				prev_mt_id = mt_id
				mt_id_counter += 1
			counter_within_one_mt_id_block += 1
			counter += 1
			writer.writerow([mt_id, prom_id, strand, bs_disp_start, bs_disp_end, core_similarity_score, \
				matrix_similarity_score, sequence, counter_within_one_mt_id_block])
			if counter%5000 == 0 and self.report:
				sys.stderr.write("%s%s/%s"%('\x08'*20, counter, mt_id_counter))
		del reader, writer
		sys.stderr.write("Done.\n")
		return counter_mt_id_ls
	
	def run(self):
		ReformatTRANSFACOutput_csvformat = '%s.csv'%self.output_fname
		self.Transform2csvFormat(self.input_dir, ReformatTRANSFACOutput_csvformat)
		counter_mt_id_ls = self.orderReformattedTRANSFACHits(ReformatTRANSFACOutput_csvformat, self.output_fname)
		os.remove(ReformatTRANSFACOutput_csvformat)
		counter_mt_id_ls.sort()
		print counter_mt_id_ls


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:br", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = None
	input_dir = None
	output_fname = None
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
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r",):
			report = 1
		
	if input_dir and output_fname:
		instance = ReformatTRANSFACOutput(hostname, dbname, schema, input_dir, output_fname, \
			debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
