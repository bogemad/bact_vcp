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


def process_bam(reads, rgid, sample_name, rgpl, rglb, base_path,temp_dir,outdir,results_bam,sample_temp_dir,log,results_dir):
	print "%s: Removing PCR duplicates" % sample_name
	replace_dir(sample_temp_dir)
	bam = os.path.join(sample_temp_dir,"%s.raw_alignment.bam" % sample_name)
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	dedup_bam = "%s.dedup_reads.bam" % sample_name
	markd_input = "INPUT=%s" % bam
	markd_output = "OUTPUT=%s" % dedup_bam
	markd_metrics = "%s.picard_md_metrics.txt" % sample_name
	markd_metrics_arg = "METRICS_FILE=%s" % markd_metrics
	realign_intervals = "%s.target.intervals" % sample_name
	preprocessed_reads = "%s.preprocessed_reads.bam" % sample_name
	bbiinput = "INPUT=%s" % dedup_bam
	cov_stats_name = "%s.coverage_stats" % sample_name
	with cd(sample_temp_dir):
		shutil.copyfile(results_bam,bam)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",markd_input],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"MarkDuplicates",markd_input,markd_output,markd_metrics_arg],log)
		print "%s: PCR duplicates removed" % sample_name
		print "%s: Realigning indels" % sample_name
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","RealignerTargetCreator","-R","../reference.fa","-I",dedup_bam,"-o",realign_intervals],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","IndelRealigner","-R","../reference.fa","-I",dedup_bam,"-targetIntervals",realign_intervals,"-o",preprocessed_reads],log)
		print "%s: Indel realignment complete." % sample_name
		print "%s: Generating depth of coverage analysis" % sample_name
		cov_stats_dir = "coverage_statistics"
		os.mkdir(cov_stats_dir)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","DepthOfCoverage","-R","../reference.fa","-o",os.path.join(cov_stats_dir,cov_stats_name),"-I",preprocessed_reads,"-ct","5","-ct","10"],log)
		print "%s: Depth of coverage analysis complete. Saved to: %s" % (sample_name, os.path.join(results_dir,cov_stats_dir))
		print "%s: Processing of bam alignment complete. Saved to: %s" % (sample_name, os.path.join(results_dir,preprocessed_reads))
	shutil.copytree(os.path.join(sample_temp_dir,cov_stats_dir),os.path.join(results_dir,cov_stats_dir))
	shutil.copy(os.path.join(sample_temp_dir,markd_metrics),os.path.join(results_dir))
	shutil.copy(os.path.join(sample_temp_dir,preprocessed_reads),os.path.join(results_dir))
	shutil.rmtree(sample_temp_dir)


if __name__ == '__main__':
	process_bam(reads, rgid, sample_name, rgpl, rglb, base_path,temp_dir,outdir,results_bam,sample_temp_dir,log,results_dir)

