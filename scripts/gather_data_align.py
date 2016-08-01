import sys, os, shutil, subprocess
from contextlib import contextmanager

base_path = os.path.dirname(os.path.realpath(__file__))
sample_key = os.path.join(base_path,"keys",os.path.basename(sys.argv[1]))
sample_name, ext = os.path.splitext(os.path.basename(sample_key))

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


def get_paired_unpaired(sample_key):
	with open(sample_key) as infile:
		pup_operator = infile.next().strip()
		if pup_operator == "paired":
			paired = True
			return paired
		elif pup_operator == "unpaired":
			paired = False
			return paired
	print "Error: %s has no pup_operator. Can't work out if reads are paired or unpaired."
	sys.exit(1)


def get_paired_dirs_and_paths(base_path,sample_key):
	with open(sample_key) as infile:
		pup_operator = infile.next().strip()
		trimmed_paired_reads1 = infile.next().strip()
		trimmed_paired_reads2 = infile.next().strip()
		trimmed_unpaired_reads1 = infile.next().strip()
		trimmed_unpaired_reads2 = infile.next().strip()
	return trimmed_paired_reads1, trimmed_paired_reads2, trimmed_unpaired_reads1, trimmed_unpaired_reads2


def get_unpaired_dirs_and_paths(base_path,sample_key):
	with open(sample_key) as infile:
		pup_operator = infile.next().strip()
		trimmed_unpaired_reads = infile.next().strip()
	return trimmed_unpaired_reads


def make_sample_outdir(outdir,sample_name):
	if os.path.isdir(os.path.join(outdir,sample_name)) == True:
		print "%s directory already exists. Two samples may have identical names. Check your sample_names and retry. Exiting..." % sample_name
		sys.exit(1)
	os.mkdir(os.path.join(outdir,sample_name))


def sam_to_bam(sam,paired,sample_name):
	script_path = os.path.dirname(os.path.realpath(__file__))
	picard_path = os.path.join(script_path,"picard-tools/picard.jar")
	bam = "%s.raw_alignment.bam" % sample_name
	bam_arg = "O=%s" % bam
	if paired == True:
		samfile_args = ["I=%s" % s for s in sam]
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"MergeSamFiles",samfile_args[0],samfile_args[1],samfile_args[2],bam_arg,"SORT_ORDER=coordinate","MERGE_SEQUENCE_DICTIONARIES=true"])
		for s in sam:
			os.remove(s)
	else:
		samfile_arg = "I=%s" % sam
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"SortSam",samfile_arg,bam_arg,"SORT_ORDER=coordinate"])
		os.remove(sam)


#read-groups!
def run_smalt(reads,paired,base_path,outdir,temp_dir,sample_name):
	smalt_path = os.path.join(base_path,"smalt/src/smalt")
	with cd(temp_dir):
		run_process([smalt_path,"index","-k","13","-s","8","reference","reference.fa"])
	if paired == True:
		paired_sam = "%s_paired.sam" % sample_name
		unpaired_sam1 = "%s_unpaired1.sam" % sample_name
		unpaired_sam2 = "%s_unpaired2.sam" % sample_name
		with cd(temp_dir):
			run_process([smalt_path,"map","-i","1500","-n","1","-o",paired_sam,"reference",reads[0],reads[1]])
			run_process([smalt_path,"map","-i","1500","-n","1","-o",unpaired_sam1,"reference",reads[2]])
			run_process([smalt_path,"map","-i","1500","-n","1","-o",unpaired_sam2,"reference",reads[3]])
			sam = paired_sam, unpaired_sam1, unpaired_sam2
	else:
		sam = "%s_unpaired.sam" % sample_name
		with cd(temp_dir):
			run_process([smalt_path,"map","-i","1500","-n","1","-o",sam,"reference",reads])
	with cd(temp_dir):
		bam = sam_to_bam(sam,paired,sample_name)
		shutil.move(bam,os.path.join(outdir,sample_name,bam))
		with open("sample.list",'a') as sample_list:
			sample_list.write("%s/n" % sample_name)


outdir = os.path.join(base_path,"results")
temp_dir = os.path.join(base_path,"intermediate_files")
make_sample_outdir(outdir,sample_name)
paired = get_paired_unpaired(sample_key)
if paired == True:
	read_files = get_paired_dirs_and_paths(base_path,sample_key)
else:
	trimmed_unpaired_reads = get_unpaired_dirs_and_paths(base_path,sample_key)
run_smalt(reads,paired,base_path,outdir,temp_dir,sample_name)


