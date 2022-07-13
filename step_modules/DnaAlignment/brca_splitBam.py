#!/usr/bin/env python2.7
import os, sys
import argparse

try: import pysam
except: sys.exit('pysam module not found.\nPlease install it before.')
try: import numpy as np
except: sys.exit('numpy module not found.\nPlease install it before.')
try: from scipy.stats import binom
except: sys.exit('scipy module not found.\nPlease install it before.')

from collections import defaultdict
from collections import OrderedDict

def get_softclip_length(read):
	if not read.qual:
		qual = np.array([40] * read.rlen)
	else:
		qual = np.array(map(ord, list(read.qual)))
		qual = qual - 33
	if read.cigar[0][0] == 4:
		if read.cigar[-1][0] == 4:
			if read.cigar[0][1] > read.cigar[-1][1]:
				return read.cigar[0][1], qual[:read.cigar[0][1]], read.pos
			else:
				return read.cigar[-1][1], qual[read.rlen-read.cigar[-1][1]:], read.aend
		else:
			return read.cigar[0][1], qual[:read.cigar[0][1]], read.pos
	elif read.cigar[-1][0] == 4:
		return read.cigar[-1][1], qual[read.rlen - read.cigar[-1][1]:], read.aend
	else:
		return 0, qual, -1

def run(args, bed, bam, exbam, rebam):
	# Get BED Files to BED list
	with open(bed, 'r') as f:
		bed_list = list()
		for line in f:
			bed_dict = OrderedDict()
			bed = line.strip().split()
			bed_dict['chr'] = bed[0]
			bed_dict['start'] = int(bed[1])
			bed_dict['end'] = int(bed[2])	
			bed_list.append(bed_dict)

	# Open Input/Output BAM file
	bam_file = pysam.Samfile(bam,'rb')
	ex_file = pysam.Samfile(exbam, 'wb', template=bam_file)
	re_file = pysam.Samfile(rebam, 'wb', template=bam_file)

	# Select and Write Read to New File (in bed file range, add original sam line number)
	extract_count = 0
	remain_count = 0
	match_flag = False

	# Parameter setting
	scliplen_cutoff = 0.2
	lowqual_cutoff = 20
	min_percent_hq = 0.8
	mapq_cutoff = 1

	total_read_cnt = int(pysam.Samfile(bam,'rb').count())

	# Iterative Search in BAM
	for i, read in enumerate(bam_file):
		match_flag = False

		# Skip unmapped & secondary alignment & duplicated read
		if read.is_secondary or read.is_unmapped or read.is_duplicate:
			match_flag = False
		else:
			match_flag = True

		# Check BED Region
		if match_flag:
			read_start = int(read.reference_start)
			read_end = int(read.reference_end)
			read_chr = read.reference_name			
			# Match with BED List
			for bed in bed_list:
				if bed['chr'] == read_chr and ((read_start >= bed['start'] and read_start <= bed['end']) or (read_end >= bed['start'] and read_end <= bed['end'])):
					match_flag = True
					break
				else:
					match_flag = False

		# Check SOFTCLIP Ratio
                '''
		if match_flag:
			soft_len, soft_qual, soft_pos = get_softclip_length(read)
			sclip_ratio = soft_len / float(read.rlen)
			if soft_pos != -1:
				sclip_hq_ratio = len(soft_qual[soft_qual >= lowqual_cutoff]) / float(len(soft_qual))
			else:
				sclip_hq_ratio = 0
			if sclip_ratio >= scliplen_cutoff and sclip_hq_ratio >= min_percent_hq and read.mapq >= mapq_cutoff:
				match_flag = True
			else:
				match_flag = False
		'''
		# Write new BAM file
		if match_flag:
			ex_file.write(read)
			extract_count += 1
		else:
			re_file.write(read)
			remain_count += 1

		# STDOUT Progress
		if args.verbose:
			progress = (float(i+1) / float(total_read_cnt) * 100)

			sys.stderr.write("[{0}% completed] {1} reads processed. {2} reads are extracted, {3} reads are remained.\r".format(int(progress), i+1, extract_count, remain_count))
			sys.stderr.flush()

	bam_file.close()
	ex_file.close()
	re_file.close()
	sys.stderr.write('\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Split BAM file (for GMAP align')
	parser.add_argument("bed", metavar="<Target BED file>", help="Target BED file")
	parser.add_argument("bam", metavar="<Input BAM file>", help="Input BAM file")
	parser.add_argument("exbam", metavar="<Output Extracted BAM file>", help="Output Extracted BAM file")
	parser.add_argument("rebam", metavar="<Output Remain BAM file>", help="Output Remain BAM file")
	parser.add_argument('-v',"--verbose", help="Verbose option", action="store_true")
	args = parser.parse_args()

	run(args, args.bed, args.bam, args.exbam, args.rebam)
