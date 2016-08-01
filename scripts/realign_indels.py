import os, subprocess
from contextlib import contextmanager

base_path = os.path.dirname(os.path.realpath(__file__))
sample_key = os.path.join(base_path,"keys",os.path.basename(sys.argv[1]))
sample_name, ext = os.path.splitext(os.path.basename(sample_key))
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
	print(output)


def realign_indels(base_path,sample_name,temp_dir,outdir):
	picard_path = os.path.join(base_path,"picard-tools/picard.jar")
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	input_bam = "%s.dedup_reads.bam" % sample_name
	realign_intervals = "%s.target_intervals.txt" % sample_name
	preprocessed_reads = "%s.gatk_preprocessed_reads.bam" % sample_name
	bbiinput = "INPUT=%s" % preprocessed_reads
	with cd(temp_dir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","RealignerTargetCreator","-R","reference.fa","-I",input_bam,"-o",realign_intervals])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","IndelRealigner","-R","reference.fa","-I",input_bam,"-targetIntervals",realign_intervals,"-o",preprocessed_reads])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput])


realign_indels(base_path,sample_name,temp_dir,outdir)
