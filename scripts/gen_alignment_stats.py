import os, shutil, subprocess
from contextlib import contextmanager

base_path = os.path.dirname(os.path.realpath(__file__))
sample_key = os.path.join(base_path,"keys",os.path.basename(sys.argv[1]))
outdir = os.path.join(base_path,"results")
temp_dir = os.path.join(base_path,"intermediate_files")
gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
sample_name = sample_key[:-4]
sample_dir = os.path.join(outdir,sample_name)
coverage_dir = os.path.join(sample_dir,"coverage_data")
if os.path.isdir(coverage_dir) == True:
	shutil.rmtree(coverage_dir)
os.mkdir(coverage_dir)
bam = "%s.gatk_preprocessed_reads.bam\n" % sample_name


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

bam_list = os.path.join(temp_dir,"%s.bamlist" % sample_name)
with open(bam_list,'w') as outfile:
	outfile.write(bam)
with cd(coverage_dir):
	run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","DepthOfCoverage","-R",os.path.join(temp_dir,"reference.fa"),"-o",sample_name,"-I",bam_list,"-ct","5","-ct","10"])