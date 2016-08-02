#!/usr/bin/env python

import os, subprocess
from contextlib import contextmanager

base_path = os.path.dirname(os.path.dirname(sys.argv[1]))
sample_key = sys.argv[1]
outdir = os.path.join(base_path,"results")
temp_dir = os.path.join(base_path,"intermediate_files")
trimo_path = os.path.join(base_path,"bin/Trimmomatic-0.36/trimmomatic-0.36.jar")
cut_adapt_path = os.path.join(base_path, ".venv/bin/cutadapt")

@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)


def run_process(cmd):
	output = subprocess.check_output(cmd, stderr = subprocess.STDOUT)
	print output


def run_piped_shell_process(cmd):
	subprocess.call(cmd, shell=True)


def interleaved(base_path,sample_key,temp_dir,trimo_path):
	with open(sample_key) as infile:
		sample_name = infile.next().strip()
		reads = infile.next().strip()
		paired_output1 = infile.next().strip()
		paired_output2 = infile.next().strip()
		unpaired_output1 = infile.next().strip()
		unpaired_output2 = infile.next().strip()
	fastqutils_path = os.path.join(base_path,"bin/ngsutils-ngsutils-0.5.9/bin/fastqutils")
	print "Importing %s: Interleaved paired reads" % sample_name 
	trim_log = "%s.trimlog" % sample_name
	ca_trim_log = "%s.ca_trimlog" % sample_name
	reads1 = "%s.1.fastq" % sample_name
	reads2 = "%s.2.fastq" % sample_name
	with cd(temp_dir):
		print "Spliting %s reads file" % sample_name 
		run_process([fastqutils_path,"unmerge",reads,sample_name])
		print "Trimming %s reads" % sample_name 
		cut_adapt_out1 = "ca_%s" % reads1
		cut_adapt_out2 = "ca_%s" % reads2
		run_piped_shell_process("%s -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o %s -p %s %s %s > %s" % (cut_adapt_path, cut_adapt_out1, cut_adapt_out2, reads1, read2, ca_trim_log))
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",trimo_path,"PE","-phred33","-trimlog",trim_log,cut_adapt_out1,cut_adapt_out2,paired_output1,unpaired_output1,paired_output2,unpaired_output2,"LEADING:20","TRAILING:20","SLIDINGWINDOW:4:15","MINLEN:70"])
		os.remove(reads1)
		os.remove(reads2)
		os.remove(cut_adapt_out1)
		os.remove(cut_adapt_out2)


interleaved(base_path,sample_name,temp_dir,outdir)
