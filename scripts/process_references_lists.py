#!/usr/bin/env python

import os, sys, shutil, subprocess
from contextlib import contextmanager


base_path = os.path.dirname(os.path.dirname(sys.argv[1]))


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


def setup_directories_and_inputs(base_path):
	outdir = os.path.join(base_path,"results")
	temp_dir = os.path.join(base_path,"intermediate_files")
	file_list = sys.argv[1]
	return outdir, temp_dir, file_list


def import_file_list(file_list):
	import_tuple = (False,False,False)
	reads_list = []
	ref_fasta = []
	ref_db = []
	with open(file_list) as file_handle:
		for line in file_handle:
			if line.startswith("#") or line == "\n":
				continue
			elif line == "[read_data]\n":
				print "Importing paired interleaved files..."
				import_tuple = (True,False,False)
				continue
			elif line == "[reference_fasta]\n":
				print "Importing reference fasta..."
				import_tuple = (False,True,False)
				continue
			elif line == "[snpEff_reference_database]\n":
				print "Importing snpEff reference database..."
				import_tuple = (False,False,True)
				continue
			elif import_tuple == (True,False,False):
				reads_list.append(line.strip())
			elif import_tuple == (False,True,False):
				ref_fasta = (line.split(',')[0])
			elif import_tuple == (False,False,True):
				ref_db = (line.split(',')[0])
	return reads_list, ref_fasta, ref_db


def copy_reference(base_path,ref_fasta,temp_dir,outdir):
	picard_path = os.path.join(base_path,"picard-tools/picard.jar")
	print "Importing reference sequence: %s" % os.path.basename(ref_fasta)
	ref_path = os.path.join(temp_dir,"reference.fa")
	shutil.copyfile(ref_fasta,os.path.join(temp_dir,"reference.fa"))
	with cd(temp_dir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"CreateSequenceDictionary","REFERENCE=reference.fa","OUTPUT=reference.dict"])
	shutil.copyfile(os.path.join(temp_dir,"reference.fa"),os.path.join(outdir,"reference.fa"))


def copy_ref_db_path(ref_db,temp_dir):
	with open(os.path.join(temp_dir,"snpEff_ref_db.txt"),'w') as outfile:
		outfile.write(ref_db)


def output_keylists(reads_list):
	for reads in reads_list:
		sample_name, file_ext = os.path.splitext(os.path.basename(reads))
		sample_outfile = "%s.txt" % sample_name
		paired_output1 = "%s.1.paired_trimmed.fastq" % sample_name
		paired_output2 = "%s.2.paired_trimmed.fastq" % sample_name
		unpaired_output1 = "%s.1.unpaired_trimmed.fastq" % sample_name
		unpaired_output2 = "%s.2.unpaired_trimmed.fastq" % sample_name
		with open(sys.argv[2],'w') as outfile:
			outfile.write("%s\n%s\n%s\n%s\n%s\n%s\n" % (sample_name,reads,paired_output1,paired_output2,unpaired_output1,unpaired_output2))



outdir, temp_dir, file_list = setup_directories_and_inputs(base_path)
reads_list, ref_fasta, ref_db = import_file_list(file_list)
null_output = copy_reference(base_path,ref_fasta,temp_dir,outdir)
null_output = copy_ref_db_path(ref_db,temp_dir)
output_keylists(reads_list)

