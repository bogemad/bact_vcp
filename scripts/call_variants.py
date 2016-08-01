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
	print output


def run_shell_process(cmd):
	output = subprocess.check_output(cmd, stderr = subprocess.STDOUT, shell=True)
	print output


def call_variants(base_path,sample_name,temp_dir,outdir):
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	preprocessed_reads = "%s.gatk_preprocessed_reads.bam" % sample_name
	bbiinput2 = "INPUT=%s" % preprocessed_reads
	raw_variants = "%s.raw_variants.vcf" % sample_name
	raw_snps = "%s.raw_snps.vcf" % sample_name
	raw_indels = "%s.raw_indels.vcf" % sample_name
	filtered_snps = "%s.filtered_snps.vcf" % sample_name
	filtered_indels = "%s.filtered_indels.vcf" % sample_name
	with cd(temp_dir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","HaplotypeCaller","-R","reference.fa","-I",preprocessed_reads,"--genotyping_mode","DISCOVERY","-glm","BOTH","-stand_call_conf","30","-stand_emit_conf","10","-o",raw_variants,"--sample_ploidy","1"])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",raw_variants,"-selectType","SNP","-o",raw_snps])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",raw_variants,"-selectType","SNP","-o",raw_indels])
		run_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s --filterExpression "QD < 3.0" --filterName "low_qual_by_depth" --filterExpression "FS > 13.5" --filterName "strand_bias" --filterExpression "DP < 10" --filterName "low_depth" --filterExpression "MQRankSum < -12.5" --filterName "MQRankSum" --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum" -o %s' % (gatk_path,raw_snps,filtered_snps))
		run_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s --filterExpression "QD < 3.0" --filterName "low_qual_by_depth" --filterExpression "FS > 200.0" --filterName "strand_bias" --filterExpression "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum" --filterExpression "DP < 10" --filterName "low_depth" -o %s' % (gatk_path,raw_indels,filtered_indels))
		shutil.copy(filtered_snps,os.path.join(outdir,sample_name))
		shutil.copy(filtered_indels,os.path.join(outdir,sample_name))


call_variants(base_path,sample_name,temp_dir,outdir)
