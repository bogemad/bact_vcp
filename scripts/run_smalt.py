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


def run_po_command(command, outputfile, errorfile):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = stdout, stderr = process.communicate()
	if output:
		with open(outputfile,'w') as outhandle:
			outhandle.write(output)
	if error:
		with open(errorfile,'a+') as errorhandle:
			errorhandle.write(error)
	if process.returncode == 1:
		print "\n\n###########    Error    ###########\n\nThe following command failed to complete successfully:\n\n%s\n\nOutput of this error can be found in:\n\n%s\n\n" % (" ".join(command), errorfile)
		sys.exit(1)


def sam_to_bam(sam, reads, rgid, sample_name, rgpl, rglb ,base_path, log):
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	bam = "%s.raw_alignment.bam" % sample_name
	pre_bam = "%s.raw_no_rg.bam" % sample_name
	run_po_command(["samtools", "view", "-bh", "-q", "10", sam], pre_bam, log)
	run_command(["java",
	"-Xmx2g","-XX:+UseSerialGC",
	"-jar",picard_path,
	"AddOrReplaceReadGroups",
	"I=%s" % pre_bam,
	"O=%s" % bam,
	"RGID=%s" % rgid,
	"RGSM=%s" % sample_name,
	"RGPL=%s" % rgpl.lower(),
	"RGLB=%s" % rglb,
	"RGPU=NA",
	"SORT_ORDER=coordinate"],log)
	os.remove(sam)
	os.remove(pre_bam)
	alignment_quality_check(bam,picard_path,sample_name,log)
	return bam


def alignment_quality_check(bam,picard_path,sample_name,log):
	stats_dir = "alignment_quality_stats"
	os.mkdir(stats_dir)
	prefix = os.path.join(os.path.abspath(stats_dir),sample_name)
	print "%s: Generating alignment quality plots" % sample_name
	run_command(["java",
	"-Xmx2g","-XX:+UseSerialGC",
	"-jar",picard_path,
	"CollectGcBiasMetrics",
	"R=../reference.fa",
	"I=%s" % bam,
	"O=%s_GCBias.txt" % prefix,
	"CHART=%s_GCBias.pdf" % prefix,
	"S=%s_summary_metrics.txt" % prefix,
	"ASSUME_SORTED=true"],log)
	run_command(["java",
	"-Xmx2g","-XX:+UseSerialGC",
	"-jar",picard_path,
	"MeanQualityByCycle",
	"R=../reference.fa",
	"I=%s" % bam,
	"O=%s_Qcycle.txt" % prefix,
	"CHART=%s_Qcycle.pdf" % prefix],log)
	run_command(["java",
	"-Xmx2g","-XX:+UseSerialGC",
	"-jar",picard_path,
	"QualityScoreDistribution",
	"R=../reference.fa",
	"I=%s" % bam,
	"O=%s_Qdist.txt" % prefix,
	"CHART=%s_Qdist.pdf" % prefix],log)
	print "%s: Alignment quality plots complete. Saved to: %s" % (sample_name, os.path.abspath(stats_dir))


def run_smalt(reads,rgid,sample_name,rgpl,rglb,sample_temp_dir,merged_trimmed,log,results_dir):
	replace_dir(sample_temp_dir)
	reads1 = "%s.1.fastq" % sample_name
	reads2 = "%s.2.fastq" % sample_name
	with cd(sample_temp_dir):
		print "%s: Processing reads for mapping" % sample_name
		unmerge(merged_trimmed, sample_name, gz=False)
		sam = "%s.sam" % sample_name
		print "%s: Mapping reads" % sample_name
		run_command(["smalt","map","-i","1500","-n","1","-o",sam,"../reference",reads1,reads2],log)
		print "%s: Read mapping complete" % sample_name
		print "%s: Generating raw alignment bam file" % sample_name
		bam = sam_to_bam(sam, reads, rgid, sample_name, rgpl, rglb ,base_path, log)
	shutil.copytree(os.path.join(sample_temp_dir,"alignment_quality_stats"),os.path.join(results_dir,"alignment_quality_stats"))
	shutil.copyfile(os.path.join(sample_temp_dir,bam),os.path.join(results_dir,bam))
	shutil.rmtree(sample_temp_dir)
	print "%s: Raw alignment bam file complete. Saved to: %s" % (sample_name, os.path.join(results_dir,bam))


if __name__ == '__main__':
	run_smalt(reads,rgid,sample_name,rgpl,rglb,sample_temp_dir,merged_trimmed,log,results_dir)

