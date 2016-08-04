#!/usr/bin/env python

import os, subprocess, sys
from contextlib import contextmanager

base_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
bam = sys.argv[1]
sample_name, ext = os.path.splitext(os.path.basename(bam))
outdir = os.path.join(base_path,"results")
temp_dir = os.path.join(base_path,"intermediate_files")

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


def dedup_reads(base_path, bam, sample_name, temp_dir,outdir):
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	markd_input = "INPUT=%s" % bam
	markd_output = "OUTPUT=%s" % sys.argv[2]
	markd_metrics = "METRICS_FILE=%s.picard_md_metrics.txt" % os.path.join(outdir,sample_name,sample_name)
	bbiinput = "INPUT=%s.dedup_reads.bam" % sample_name
	with cd(temp_dir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"MarkDuplicates",markd_input,markd_output,markd_metrics])
	with cd(outdir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput])


dedup_reads(base_path, bam, sample_name, temp_dir,outdir)
