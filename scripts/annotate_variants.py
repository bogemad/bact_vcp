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
	print(output)


def run_piped_shell_process(cmd):
	subprocess.call(cmd, shell=True)


def annotate_variants(base_path,sample_name,temp_dir,outdir,ref_db):
	snpeff_path = os.path.join(script_path,"snpEff/snpEff.jar")
	filtered_snps = "%s.filtered_snps.vcf" % sample_name
	filtered_indels = "%s.filtered_indels.vcf" % sample_name
	annotated_snps = "%s.annotated_snps.vcf" % sample_name
	annotated_indels = "%s.annotated_indels.vcf" % sample_name
	with cd(temp_dir):
		run_piped_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -v %s %s > %s' % (snpeff_path,ref_db,filtered_snps,annotated_snps))
		run_piped_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -v %s %s > %s' % (snpeff_path,ref_db,filtered_indels,annotated_indels))
		shutil.copy(annotated_snps,os.path.join(outdir,sample_name))
		shutil.copy(annotated_indels,os.path.join(outdir,sample_name))


annotate_variants(base_path,sample_name,temp_dir,outdir)
