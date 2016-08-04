#!/usr/bin/env python

import os, subprocess
from contextlib import contextmanager

base_path = os.path.dirname(os.path.realpath(__file__))
vcf = sys.argv[1]
sample_name, ext = os.path.splitext(os.path.basename(vcf))
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


def run_piped_shell_process(cmd):
	subprocess.call(cmd, shell=True)


def annotate_variants(base_path,sample_name,temp_dir,outdir,ref_db,vcf):
	gatk_path = os.path.join(base_path,"bin/gatk/GenomeAnalysisTK.jar")
	raw_indels = "%s.raw_indels.vcf" % sample_name
	filtered_indels = "%s.filtered_indels.vcf" % sample_name
	snpeff_path = os.path.join(script_path,"snpEff/snpEff.jar")
	annotated_indels = "%s.annotated_indels.vcf" % sample_name
	with cd(temp_dir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",vcf,"-selectType","SNP","-o",raw_indels])
		run_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s --filterExpression "QD < 3.0" --filterName "low_qual_by_depth" --filterExpression "FS > 200.0" --filterName "strand_bias" --filterExpression "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum" --filterExpression "DP < 10" --filterName "low_depth" -o %s' % (gatk_path,raw_indels,filtered_indels))
		run_piped_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -v %s %s > %s' % (snpeff_path,ref_db,filtered_indels,annotated_indels))
		shutil.copy(filtered_indels,os.path.join(outdir,sample_name))


annotate_variants(base_path,sample_name,temp_dir,outdir)
