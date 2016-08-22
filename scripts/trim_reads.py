#!/usr/bin/env python

import sys, os, shutil, subprocess
from contextlib import contextmanager


base_path = sys.argv[1]
sample_name = sys.argv[2]
reads = sys.argv[3]
rgid = sys.argv[4]
rgpl = sys.argv[5]
rglb = sys.argv[6]
results_dir = sys.argv[7]
sample_temp_dir = sys.argv[8]
log = sys.argv[9]
merged_trimmed = sys.argv[10]
bam = sys.argv[11]
dedup_bam = sys.argv[12]
preprocessed_reads = sys.argv[13]
raw_variants = sys.argv[14]


@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)


def replace_dir(dir):
	if os.path.exists(dir):
		shutil.rmtree(dir)
	os.mkdir(dir)


def run_command(command, outputfile):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	while True:
		output = process.stdout.readline()
		if output == '' and process.poll() is not None:
			break
		if output:
			with open(outputfile,'a+') as outhandle:
				outhandle.write(output)
	if process.poll() == 1:
		print "\n\n###########    Error    ###########\n\nThe following command failed to complete successfully:\n\n%s\n\nOutput of this error can be found in:\n\n%s\n\n" % (" ".join(command), outputfile)
		sys.exit(1)


def trim_reads(base_path, sample_temp_dir, sample_name, log, reads, results_dir):
	trimo_path = os.path.join(base_path,"bin/trimmomatic/trimmomatic-0.36.jar")
	replace_dir(sample_temp_dir)
	paired_output1 = "%s.trimmed_paired1.fastq" % sample_name
	paired_output2 = "%s.trimmed_paired2.fastq" % sample_name
	unpaired_output1 = "%s.trimmed_unpaired1.fastq" % sample_name
	unpaired_output2 = "%s.trimmed_unpaired2.fastq" % sample_name
	merged_trimmed = "%s.trimmed.fastq.gz" % sample_name
	print "%s: Importing reads" % sample_name 
	trim_log = "%s.trimlog" % sample_name
	reads1 = "%s.1.fastq" % sample_name
	reads2 = "%s.2.fastq" % sample_name
	with cd(sample_temp_dir):
		print "%s: Processing reads file for trimming" % sample_name 
		unmerge(reads, sample_name, gz=False)
		print "%s: Trimming reads" % sample_name 
		cut_adapt_out1 = "ca_%s" % reads1
		cut_adapt_out2 = "ca_%s" % reads2
		run_command(["cutadapt",
		"-g","AGATGTGTATAAGAGACAG",
		"-a","CTGTCTCTTATACACATCT",
		"-g","AGATGTGTATAAGAGACAG",
		"-a","CTGTCTCTTATACACATCT",
		"-g","TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
		"-a","CTGTCTCTTATACACATCTGACGCTGCCGACGA",
		"-g","GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
		"-a","CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
		"-o",cut_adapt_out1,"-p",cut_adapt_out2,
		reads1, reads2],log)
		run_command(["java",
		"-Xmx2g","-XX:+UseSerialGC",
		"-jar",trimo_path,"PE","-phred33",
		"-threads","1",
		"-trimlog",trim_log,
		cut_adapt_out1,cut_adapt_out2,
		paired_output1,unpaired_output1,
		paired_output2,unpaired_output2,
		"LEADING:20","TRAILING:20",
		"SLIDINGWINDOW:4:15","MINLEN:70"],log)
		print "%s: Merging trimmed reads" % sample_name
		with gzip.open(merged_trimmed, 'wb') as trimmed_merged:
			fastqs = [FASTQ(x) for x in (paired_output1,paired_output2)]
			merge(fastqs, split_slashes=False, out=trimmed_merged, quiet=False)
	shutil.copy(os.path.join(sample_temp_dir,merged_trimmed),os.path.join(results_dir,merged_trimmed))
	shutil.rmtree(sample_temp_dir)
	print "%s: Read trimming complete. Saved to: %s" % (sample_name, os.path.join(results_dir,merged_trimmed))


if __name__ == '__main__':
	trim_reads(base_path, sample_temp_dir, sample_name, log, reads, results_dir)
