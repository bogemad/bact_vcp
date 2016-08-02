#!/usr/bin/env python

import os, subprocess
from contextlib import contextmanager

base_path = os.path.dirname(os.path.realpath(__file__))
sample_key = os.path.join(base_path,"keys",os.path.basename(sys.argv[1]))
sample_name, ext = os.path.splitext(os.path.basename(sample_key))
outdir = os.path.join(base_path,"results")
temp_dir = os.path.join(base_path,"intermediate_files")
with open(os.path.join(temp_dir,"snpEff_ref_db.txt")) as infile:
	ref_db = infile.readline().strip()


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


def run_piped_shell_process(cmd):
	subprocess.call(cmd, shell=True)


def annotate_variants(base_path,sample_name,temp_dir,outdir,ref_db):
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	raw_variants = "%s.raw_variants.vcf" % sample_name
	raw_snps = "%s.raw_snps.vcf" % sample_name
	filtered_snps = "%s.filtered_snps.vcf" % sample_name
	snpeff_path = os.path.join(script_path,"snpEff/snpEff.jar")
	annotated_snps = "%s.annotated_snps.vcf" % sample_name
	with cd(temp_dir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",raw_variants,"-selectType","SNP","-o",raw_snps])
		run_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s --filterExpression "QD < 3.0" --filterName "low_qual_by_depth" --filterExpression "FS > 13.5" --filterName "strand_bias" --filterExpression "DP < 10" --filterName "low_depth" --filterExpression "MQRankSum < -12.5" --filterName "MQRankSum" --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum" -o %s' % (gatk_path,raw_snps,filtered_snps))
		run_piped_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -v %s %s > %s' % (snpeff_path,ref_db,filtered_snps,annotated_snps))
		shutil.copy(annotated_snps,os.path.join(outdir,sample_name))
		shutil.copy(annotated_indels,os.path.join(outdir,sample_name))


annotate_variants(base_path,sample_name,temp_dir,outdir)
