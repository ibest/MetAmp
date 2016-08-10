# -*- coding:utf8 -*-
#!/usr/bin/env python


"""
This script extracts variable regions from 16S rRNA. Author: Ilya Y. Zhbannikov, 2015

Example run:

python extract.py -i gold.fa -o test -fp CCTACGGGAGGCAGCAG -rp CCGTCAATTCMTTTRAGN

Here: 
gold.fa - an input file that contains 16S sequences, must be provided in FASTA format;
test - an output prefix for file that contains extracted variable regions, will be in FASTA format (i.e. test.fasta);
CCTACGGGAGGCAGCAG - forward primer
CCGTCAATTCMTTTRAGN - reverse primer
n - number of records to extract, 100 by default


Some primer sequences that can be used
V1-3 primers:
forward_primer = "NGAGTTTGATCCTGGCTCAG" # -m=2 -p=-5 g=-15
reverse_primer = "ATTACCGCGGCTGCTGG"

V3-5 primers:
forward_primer = "CCTACGGGAGGCAGCAG" # -m=2 -p=-5 g=-8
reverse_primer = "CCGTCAATTCMTTTRAGN"

V6-9 primers:
forward_primer = "GAATTGACGGGGRCCC" # -m=2 -p=-5 -g=-1
reverse_primer = "TACGGYTACCTTGTTAYGACTT"
 

Last change was made on 2014/10/14

"""

import sys
from Bio import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
#from optparse import OptionParser
import argparse
from smgb import *
from AmbiguityTable import *


class Extract :
	match = None # Match award
	mismatch = None # Mismatch penalty
	gap = None # Gap penalty
	forward_primer = None
	reverse_primer = None
	start_avg = []
	stop_avg = []
	length_avg = []
	pos_trim = False
	spos = 0
	epos = 1600
	num_reads = 3000
	verbose = False
	
	def __init__(self): # Initialization
		self.match = 2 # Matсh award (-ma)
		self.mismatch = -5 # Gap mismatch penalty (-mp)
		self.gap = -8 # (-gap) Gap penalty
		self.forward_primer = "NGAGTTTGATCCTGGCTCAG"
		self.reverse_primer = "ATTACCGCGGCTGCTGG"
		
 
	def extract(self, bacteria_filename, output_prefix) :
		# Main method that extracts variable regions from complete 16S sequences
		handle = open(bacteria_filename, "rU")
		#Output files:
		report = output_prefix + ".txt"
		
		report_handle = open(report, "w")
		not_found_handle = open(output_prefix+"_not_found.fasta", 'w')
		output_table_handle = open(output_prefix+"_table.txt",'w')
		output_prefix = output_prefix+".fasta"
		output_handle = open(output_prefix, "w")
		
		print "Output file:", output_prefix
		#print "Records to extract:", self.num_rec
		
		prev_start = 0
		prev_stop = 0
		start_pos = stop_pos = 0
		# For each line:
		cnt = 0 # Counter
		for record in SeqIO.parse(handle, "fasta") :
			
			#Try to find a forward primer
			v_region_found = False
			# We generate primers from ambiguity table:
			forward_primers_set = Primer().get_primer_set(self.forward_primer)
			start = Smgb()
			start.set_params(self.match,self.mismatch,self.gap)
			scores = []
			positions = []
			# For each generated primer (we are using an ambiguity table):
			for p in forward_primers_set :
				results = start.smgb(record.seq.upper(),p)
				positions.append(results[0])
				scores.append(results[1])
			start_pos = positions[scores.index(max(scores))]
			if start_pos != -1 : # Start found! Let's try to find stop position
				start_pos += len(self.forward_primer)
				#Try to find a reverse primer:
				reverse_primers_set = Primer().get_primer_set(self.reverse_primer)
				stop = Smgb()
				stop.set_params(self.match,self.mismatch,self.gap)
				scores = []
				positions = []
				for p in reverse_primers_set :
					results = stop.smgb(record.seq.upper(),Seq(self.reverse_primer, generic_dna).reverse_complement())
					positions.append(results[0])
					scores.append(results[1])
					stop_pos = positions[scores.index(max(scores))]
				if stop_pos > start_pos + 2*len(self.forward_primer) : 
					SeqIO.write(record[start_pos:stop_pos].upper(), output_handle, "fasta")
					v_region_found = True
					prev_start = start_pos
					prev_stop = stop_pos
			if v_region_found == False :
				forward_primers_set = Primer().get_primer_set(self.forward_primer)
				start = Smgb()
				start.set_params(self.match,self.mismatch,self.gap)
				positions = []
				scores = []
				# Try reverse complement
				revcomp = record.reverse_complement()
				for p in forward_primers_set :
					results = start.smgb(revcomp.seq.upper(),p)
					positions.append(results[0])
					scores.append(results[1])
				start_pos = positions[scores.index(max(scores))]
				if start_pos != -1 : # Start found! Let's try to find stop position
					start_pos += len(self.forward_primer)
					#Try to find a reverse primer:
					reverse_primers_set = Primer().get_primer_set(self.reverse_primer)
					stop = Smgb()
					stop.set_params(self.match,self.mismatch,self.gap)
					scores = []
					positions = []
					for p in reverse_primers_set :
						results = stop.smgb(revcomp.seq.upper(),Seq(self.reverse_primer, generic_dna).reverse_complement())
						positions.append(results[0])
						scores.append(results[1])
						stop_pos = positions[scores.index(max(scores))]
					if stop_pos > start_pos + 2*len(self.forward_primer) :
						new_record = SeqRecord(Seq(str(revcomp[start_pos:stop_pos].reverse_complement().seq.upper()), generic_dna), id=record.id, name=record.name, description=record.description)
						SeqIO.write(new_record, output_handle, "fasta")
						v_region_found = True
						prev_start = start_pos
						prev_stop = stop_pos
			if v_region_found == False : # Report the record for which marker region was not found:
				SeqIO.write(record, not_found_handle, "fasta")
			
			if self.verbose :
				print record.id, start_pos, stop_pos, v_region_found
			
			
			output_table_handle.write(record.id + "\t" + str(start_pos) + '\t' + str(stop_pos) + '\t' + str(v_region_found) + '\n')
			self.start_avg.append(start_pos)
			self.stop_avg.append(stop_pos)
			self.length_avg.append(stop_pos-start_pos)
			
			cnt += 1
			
			if cnt == self.num_reads :
				print "Done." 
				print self.num_reads, "records extracted."
				break
			
		# Cleaning up:
		handle.close()
		report_handle.write("Total records\tAvg start pos\tAvg end pos\tAvg length\n")
		report_handle.write(str(cnt)+'\t'+str(sum(self.start_avg)/len(self.start_avg))+'\t'+str(sum(self.stop_avg)/len(self.stop_avg))+'\t'+str(sum(self.length_avg)/len(self.length_avg))+'\n')
		output_handle.close()  
		report_handle.close()
		output_table_handle.close()
		
	def positional(self, filename, output_prefix) :
		handle = open(filename, "rU")
		
		#Output files:
		report = output_prefix + ".txt"
		report_handle = open(report, "w")
		output_prefix = output_prefix+".fasta"
		output_handle = open(output_prefix, "w")
		
		records = []
		i=0
		for record in SeqIO.parse(handle, "fasta") :
			#if i < self.num_reads/2 :
			record.seq = record.seq[self.spos:self.epos]
			records.append(record)
			#else :
			#	revcomp = record.reverse_complement()
			#	new_record = SeqRecord(Seq(str(revcomp[self.spos:self.epos].reverse_complement().seq.upper()), generic_dna), id=record.id, name=record.name, description=record.description)
			#	SeqIO.write(new_record, output_handle, "fasta")
			i+=1
		
		SeqIO.write(records, output_handle, "fasta")
		output_handle.close()
		handle.close()
		report_handle.write("Total records processed:\t")
		report_handle.write(str(i)+'\n')
		output_handle.close()  
		report_handle.close()
		
    	
def main() :
	# Setting-up options for option parser:
	usage = "extract.py -i <input file path> -o <output file prefix> [options]"
	parser = argparse.ArgumentParser(usage=usage)
	parser.add_argument("-i", "--input", action="store", type=str, dest="infile", default=None, required=True, help="input file") # V13 forward primer
	parser.add_argument("-o", "--output", action="store", type=str, dest="outprefix", default=None, required=True, help="output prefix") # V13 forward primer
	
	parser.add_argument("-fp", "--forward_primer", action="store", type=str, dest="forward_primer", default=None, help="forward primer sequence") # V13 forward primer
	parser.add_argument("-rp", "--reverse_primer", action="store", type=str, dest="reverse_primer", default=None, help="reverse primer sequence") # V13 reverse primer
	parser.add_argument("-m", "--match_award", action="store", type=int, dest="match_award", default=2, help="match award")
	parser.add_argument("-p", "--mismatch_penalty", action="store", type=int, dest="mismatch_penalty", default=-5, help="mismatch penalty")
	parser.add_argument("-g", "--gap", action="store", type=int, dest="gap_penalty", default=-8, help="gap penalty")
	## Positional trimming: the user has to set the coordinates manually.
	parser.add_argument("-pos", "--positional", action="store_true", dest="pos_trim", default=False, help="positional trimming")
	parser.add_argument("-spos", "--start_position", action="store", type=int, dest="spos", default=0, help="positional trimming: from")
	parser.add_argument("-epos", "--end_position", action="store", type=int, dest="epos", default=1600, help="positional trimming: to")
	parser.add_argument("-nr", "--num_reads", action="store", type=int, dest="num_reads", default=3000, help="number of reads in a library")
	parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False, help="verbosing output")
	
	args = parser.parse_args()
	
	# Checking the input parameters:
	if len(sys.argv) < 3 :
		print "Error: no input files provided."
		parser.print_help()
	else :
		if (args.forward_primer == None or args.reverse_primer == None) and args.pos_trim == False :
			print "Warning: no primer sequences provided. The program will use V1-3 forward and reverse primer sequences: {NGAGTTTGATCCTGGCTCAG, ATTACCGCGGCTGCTGG}"
			args.forward_primer = "NGAGTTTGATCCTGGCTCAG"
			args.reverse_primer = "ATTACCGCGGCTGCTGG"
			
		# Extracter instantiation:
		obj = Extract()
		obj.match = args.match_award
		obj.mismatch = args.mismatch_penalty
		obj.gap = args.gap_penalty
		obj.forward_primer = args.forward_primer
		obj.reverse_primer = args.reverse_primer
		obj.pos_trim = args.pos_trim
		obj.spos = args.spos
		obj.epos = args.epos
		obj.num_reads = args.num_reads
		obj.verbose = args.verbose
		
		print "===========Input parameters=========="
		print "Input file:", args.infile
		print "Output prefix:", args.outprefix
		print "Match award:", args.match_award
		print "Mismatch penalty:", args.mismatch_penalty
		print "Gap:", args.gap_penalty
		print "Forward primer", args.forward_primer
		print "Reverse primer", args.reverse_primer
		
		if args.pos_trim == True :
			print "Positional trimming:", "True"
			print "Start position:", args.spos
			print "End position:", args.epos
		
		print "Reads to extract:", args.num_reads
		print "Verbose:", "True" if args.verbose == True else "False"
		print "====================================="
		
		print "Extraction started..."
		if obj.pos_trim == True :
			print "Positional trimming..."
			if obj.spos == None :
				print "Error: starting position is None. Exiting..."
				return
			elif obj.epos == None :
				print "Error: ending position is None. Exiting..."
				return
			obj.positional(args.infile, args.outprefix)
		else :
			obj.extract(args.infile, args.outprefix)
	
if __name__ == '__main__':
    main()
    
