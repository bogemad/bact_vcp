#!/usr/bin/env python

import os, subprocess
from contextlib import contextmanager

base_path = os.path.dirname(os.path.realpath(__file__))
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


def run_shell_process(cmd):
	output = subprocess.check_output(cmd, stderr = subprocess.STDOUT, shell=True)
	print output


def call_variants(base_path,sample_name,temp_dir,outdir,bam):
	gatk_path = os.path.join(base_path,"bin/gatk/GenomeAnalysisTK.jar")
	raw_variants = sys.argv[2]
	with cd(temp_dir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","HaplotypeCaller","-R","reference.fa","-I",bam,"--genotyping_mode","DISCOVERY","-glm","BOTH","-stand_call_conf","30","-stand_emit_conf","10","-o",raw_variants,"--sample_ploidy","1"])


call_variants(base_path,sample_name,temp_dir,outdir,bam)
