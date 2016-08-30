#!/usr/bin/env python

import os
from multiprocessing import Pool
from bvcp_core import trim_reads, run_smalt, process_bam, call_variants

def run_pipeline_local(threads, reads_list):
	# for read_data in reads_list:
		# multiprocessing(read_data)
	p = Pool(int(threads))
	p.map_async(multiprocessing,reads_list).get(9999999)
	p.close()
	p.join()


def run_pipeline(read_data,base_path,temp_dir):
	outdir = read_data[5]
	results_dir = os.path.join(outdir,read_data[2])
	if os.path.exists(results_dir) == False:
		os.mkdir(results_dir)
	sample_temp_dir = os.path.join(temp_dir,read_data[2])
	log = "%s/%s.log" % (results_dir,read_data[2])
	merged_trimmed = os.path.join(results_dir,"%s.trimmed.fastq.gz" % read_data[2])
	bam = os.path.join(results_dir,"%s.raw_alignment.bam" % read_data[2])
	dedup_bam = os.path.join(results_dir,"%s.dedup_reads.bam" % read_data[2])
	preprocessed_reads = os.path.join(results_dir,"%s.preprocessed_reads.bam" % read_data[2])
	raw_variants = os.path.join(results_dir,"%s.raw_variants.g.vcf" % read_data[2])
	if os.path.exists(merged_trimmed) == False:
		trim_reads(read_data,base_path,temp_dir,sample_temp_dir,log,results_dir)
	if os.path.exists(bam) == False:
		run_smalt(read_data,base_path,temp_dir,merged_trimmed,sample_temp_dir,log,results_dir)
	if os.path.exists(preprocessed_reads) == False:
		process_bam(read_data,base_path,temp_dir,bam,sample_temp_dir,log,results_dir)
	if os.path.exists(raw_variants) == False:
		call_variants(read_data,base_path,temp_dir,preprocessed_reads,sample_temp_dir,log,results_dir)


def multiprocessing(read_data):
	base_path = os.path.dirname(os.path.dirname(__file__))
	temp_dir = os.path.join(base_path,"intermediate_files")
	run_pipeline(read_data,base_path,temp_dir)

