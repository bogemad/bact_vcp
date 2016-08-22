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
results_bam = sys.argv[11]
dedup_bam = sys.argv[12]
results_preprocessed_reads = sys.argv[13]
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

def call_variants(reads, rgid, sample_name, rgpl, rglb,base_path,temp_dir,outdir,results_preprocessed_reads,sample_temp_dir,log,results_dir):
	print "%s: Calling variants" % sample_name
	replace_dir(sample_temp_dir)
	preprocessed_reads = "%s.preprocessed_reads.bam" % sample_name
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	raw_variants = "%s.raw_variants.g.vcf" % sample_name
	bbiinput = "INPUT=%s" % preprocessed_reads
	with cd(sample_temp_dir):
		shutil.copyfile(results_preprocessed_reads,preprocessed_reads)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","HaplotypeCaller","-R","../reference.fa","-I",preprocessed_reads,"--genotyping_mode","DISCOVERY","-stand_call_conf","30","-stand_emit_conf","10","-ERC", "GVCF", "-o",raw_variants,"--sample_ploidy","1"],log)
	shutil.copyfile(os.path.join(sample_temp_dir,raw_variants),os.path.join(results_dir,raw_variants))
	shutil.rmtree(sample_temp_dir)
	print "%s: Variant calling complete. Raw variants vcf saved to: %s" % (sample_name, os.path.join(results_dir,raw_variants))


if __name__ == '__main__':
	call_variants(reads, rgid, sample_name, rgpl, rglb,base_path,temp_dir,outdir,results_preprocessed_reads,sample_temp_dir,log,results_dir)

