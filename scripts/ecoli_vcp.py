#Runs variant calling pipeline

from contextlib import contextmanager
import os, argparse, subprocess, sys, shutil
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Runs a variant calling pipeline')

parser.add_argument('--intreads', '-p', action='append', help='One or several interleaved paired read sequencing files in .fastq format. File extension must be .fastq or .fq. Use this flag when paired reads are supplied as 1 files per sample. Use the -p flag before each file. If the reads are supplied as 2 files per sample use the -1 and -2 arguments.')
parser.add_argument('--unpaired', '-u',action='append',help="One or several unpaired sequencing read files in .fastq format. Use the -u flag before each file.")
parser.add_argument('--reference','-r', nargs=1, help='Reference genome in fasta format')
parser.add_argument('--reads1', '-1', action='append', help='One or several paired sequencing read files in .fastq format. Use this flag when paired reads are supplied as 2 files per sample. Use the -1 flag before each file that is first of a pair. Filenames must end in 1.fastq or 1.fq with the remaining filename identical to the other paired read file.')
parser.add_argument('--reads2', '-2', action='append', help='One or several paired sequencing read files in .fastq format. Use this flag when paired reads are supplied as 2 files per sample. Use the -2 flag before each file that is first of a pair. Filenames must end in 2.fastq or 2.fq with the remaining filename identical to the other paired read file.')
parser.add_argument('--threads', '-t', nargs=1, type=int, help='Number of threads for multiprocessing.')
parser.add_argument('--output_dir', '-o', nargs=1, help='Directory for output files')
parser.add_argument('--keep_temporary_files', '-k', action='store_true', help='Flag to keep temporary files, if this argument is included temporary files will included with regular output files in each sample folder.')
parser.add_argument('--snpeff_ref_db','-d', nargs=1, help='Database name of reference genome for snpEff analysis. Find with command: java -jar <path_to_vcp_directory>/snpEff/snpEff.jar databases')



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


def run_shell_process(cmd):
	output = subprocess.check_output(cmd, stderr = subprocess.STDOUT, shell=True)
	print(output)


def run_piped_shell_process(cmd):
	subprocess.call(cmd, shell=True)



def make_out_temp_dir(arg):
	outdir = arg.output_dir
	count = 1
	while os.path.isdir(outdir) == True:
		if count == 1:
			outdir = "%s_%d" % (outdir,count)
			print("Previous output directory detected...")
		else:
			outdir = "%s_%d" (outdir[:-2],count)
		count += 1
	if count != 1:
		print("Output directory renamed to %s" % outdir)
	os.mkdir(outdir)
	temp_dir = os.path.join(outdir,"tmp")
	os.mkdir(temp_dir)
	return outdir, temp_dir


def match_reads1_reads2(arg,read1_filename,read1_path):
	for read2_path in arg.reads2:
		read2_file = os.path.basename(read2_path)
		read2_filename, read2_file_extension = os.path.splitext(read2_file)
		if read1_filename[:-1] == read2_filename[:-1]:
			return read1_path, read2_path
	print("%s has no matching paired reads file" % read1_path)
	sys.exit(1)


def split_interleaved(temp_dir, arg):
	readpath_d = {}
	for file in arg.intreads:
		filename, file_extension = os.path.splitext(file)
		if file_extension != .fastq or file_extension != .fq:
			print("%s is not a fastq file. Ensure that all fastq files are correctly labeled with a .fq or .fastq file extension." % os.path.basename(file))
		sample_name = os.path.basename(filename)
		with cd(temp_dir):
			run_process(["fastqutils","unmerge",file,sample_name])
			fastq1 = "%s.1.fastq" % sample_name
			fastq2 = "%s.2.fastq" % sample_name
		readpath_d[sample_name] = (os.path.join(temp_dir,fastq1),os.path.join(temp_dir,fastq2))
	return readpath_d


def gen_basename_indexed_read_path_dict(arg):
	readpath_d = split_interleaved(temp_dir, arg)
	for read1_path in arg.reads1:
		read1_file = os.path.basename(read1_path)
		read1_filename, read1_file_extension = os.path.splitext(read1_file)
		sample_name = read1_filename[:-1].rstrip(["._-")
		readpath_d[sample_name] = match_reads1_reads2(arg,read1_filename,read1_path)
	return readpath_d


def sam_to_bam(sam):
	script_path = os.path.dirname(os.path.realpath(__file__))
	samtools_path = os.path.join(script_path,"samtools/samtools")
	bam = "%s.raw_alignment.bam" % sam[:-4]
	run_piped_shell_process(["%s view -bh -q 10 %s > %s" % (samtools_path,sam,bam)])
	os.remove(sam)
	return bam


def make_sample_outdir(outdir,sample_name):
	if os.path.isdir(os.path.join(outdir,sample_name)) == True:
		print("%s directory already exists. Two samples may have identical names. Check your sample_names and retry. Exiting..." % sample_name)
		sys.exit(1)
	os.mkdir(os.path.join(outdir,sample_name))


def copy_ref(arg,temp_dir,outdir):
	saved_ref_path = os.path.join(outdir,os.path.basename(arg.reference[0]))
	temp_ref_path = os.path.join(temp_dir,"reference.fa")
	shutil.copyfile(arg.reference[0],saved_ref_path)
	shutil.copyfile(arg.reference[0],temp_ref_path)


def run_smalt(read1_path,read2_path,outdir,temp_dir,sample_name):
	script_path = os.path.dirname(os.path.realpath(__file__))
	smalt_path = os.path.join(script_path,"smalt/src/smalt")
	with cd(temp_dir):
		run_process([smalt_path,"index","-k","13","-s","8","reference","reference.fa"])
	sam = "%s.sam" % sample_name
	with cd(temp_dir):
		run_process([smalt_path,"map","-i","1500","-n","1","reference",read1_path,read2_path,sam])
		bam = sam_to_bam(sam)
		shutil.move(bam,os.path.join(outdir,sample_name,bam))


def call_variants(sample_name,temp_dir,outdir,ref_db):
	script_path = os.path.dirname(os.path.realpath(__file__))
	picard_path = os.path.join(script_path,"picard-tools/picard.jar")
	gatk_path = os.path.join(script_path,"gatk/GenomeAnalysisTK.jar")
	snpeff_path = os.path.join(script_path,"snpEff/snpEff.jar")
	bam = "%s.raw_alignment.bam" % sample_name
	sortsaminput = "INPUT=%s" % os.path.join(outdir,sample_name,bam)
	sortsamoutput = "OUTPUT=%s.sorted_reads.bam" % sample_name
	markd_input = "INPUT=%s.sorted_reads.bam" % sample_name
	markd_output = "OUTPUT=%s.dedup_reads.bam" % sample_name
	markd_metrics = "METRICS_FILE=%s.picard_md_metrics.txt" % os.path.join(outdir,sample_name,sample_name)
	bbiinput1 = "INPUT=%s.dedup_reads.bam" % sample_name
	input_bam = "%s.dedup_reads.bam" % sample_name
	realign_intervals = "%s.target_intervals.txt" % sample_name
	preprocessed_reads = "%s.gatk_preprocessed_reads.bam" % sample_name
	bbiinput2 = "INPUT=%s" % preprocessed_reads
	raw_variants = "%s.raw_variants.vcf" % sample_name
	raw_snps = "%s.raw_snps.vcf" % sample_name
	raw_indels = "%s.raw_indels.vcf" % sample_name
	filtered_snps = "%s.filtered_snps.vcf" % sample_name
	filtered_indels = "%s.filtered_indels.vcf" % sample_name
	annotated_snps = "%s.annotated_snps.vcf" % sample_name
	annotated_indels = "%s.annotated_indels.vcf" % sample_name
	with cd(temp_dir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"CreateSequenceDictionary","REFERENCE=reference.fa","OUTPUT=reference.dict"])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"SortSam",sortsaminput,sortsamoutput,"SORT_ORDER=coordinate"])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"MarkDuplicates",markd_input,markd_output,markd_metrics])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput1])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","RealignerTargetCreator","-R","reference.fa","-I",input_bam,"-o",realign_intervals])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","IndelRealigner","-R","reference.fa","-I",input_bam,"-targetIntervals",realign_intervals,"-o",preprocessed_reads])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput2])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","HaplotypeCaller","-R","reference.fa","-I",preprocessed_reads,"--genotyping_mode","DISCOVERY","-glm","BOTH","-stand_call_conf","30","-stand_emit_conf","10","-o",raw_variants,"--sample_ploidy","1"])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",raw_variants,"-selectType","SNP","-o",raw_snps])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",raw_variants,"-selectType","SNP","-o",raw_indels])
		run_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s --filterExpression "QD < 3.0" --filterName "low_qual_by_depth" --filterExpression "FS > 13.5" --filterName "strand_bias" --filterExpression "DP < 10" --filterName "low_depth" --filterExpression "MQRankSum < -12.5" --filterName "MQRankSum" --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum" -o %s' % (gatk_path,raw_snps,filtered_snps))
		run_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s --filterExpression "QD < 3.0" --filterName "low_qual_by_depth" --filterExpression "FS > 200.0" --filterName "strand_bias" --filterExpression "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum" --filterExpression "DP < 10" --filterName "low_depth" -o %s' % (gatk_path,raw_indels,filtered_indels))
		run_piped_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -v %s %s > %s' % (snpeff_path,ref_db,filtered_snps,annotated_snps))
		run_piped_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -v %s %s > %s' % (snpeff_path,ref_db,filtered_indels,annotated_indels))
		shutil.copy(filtered_snps,os.path.join(outdir,sample_name))
		shutil.copy(filtered_indels,os.path.join(outdir,sample_name))
		shutil.copy(annotated_snps,os.path.join(outdir,sample_name))
		shutil.copy(annotated_indels,os.path.join(outdir,sample_name))


def multiprocess_pipeline(multiprocessing_arguments):
	read1_path = multiprocessing_arguments[0]
	read2_path = multiprocessing_arguments[1]
	outdir = multiprocessing_arguments[2]
	temp_dir = multiprocessing_arguments[3]
	sample_name = multiprocessing_arguments[4]
	ref_db = multiprocessing_arguments[5]
	make_sample_outdir(outdir,sample_name)
	run_smalt(read1_path,read2_path,outdir,temp_dir,sample_name)
	# add an alignment validation step (depthofcoverage,etc)
	call_variants(sample_name,temp_dir,outdir,ref_db)



def keep_temp_files(sample_name, temp_dir, outdir):
	for filename in (f for f in os.listdir(temp_dir) if os.path.isfile(f) == True):
		if sample_name in filename:
			with cd(temp_dir):
				shutil.copy(filename,os.path.join(outdir,sample_name))


if __name__ == '__main__':
	arg = parser.parse_args()
	threads = arg.threads[0]
	outdir, temp_dir = make_out_temp_dir(arg)
	readpath_d = gen_basename_indexed_read_path_dict(arg)(temp_dir, arg)
	ref_path = copy_ref(arg,outdir)
	ref_db = arg.snpeff_ref_db[0]
	multiprocessing_arguments = [(readpath_d[s][0],readpath_d[s][1],outdir,temp_dir,s,ref_db) for s in list(readpath_d)]
	p = Pool(threads)
	p.imap(multiprocess_pipeline,multiprocessing_arguments)
	p.close()
	p.join()
	if arg.keep_temporary_files == True:
		for sample_name in list(readpath_d):
			keep_temp_files(sample_name, temp_dir, outdir)




