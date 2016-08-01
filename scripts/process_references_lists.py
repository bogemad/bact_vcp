import os, shutil, subprocess
from contextlib import contextmanager


base_path = os.path.dirname(os.path.realpath(__file__))


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


def setup_directories_and_inputs(base_path):
	outdir = os.path.join(base_path,"results")
	temp_dir = os.path.join(base_path,"intermediate_files")
	file_list = os.path.join(base_path,"input_data_list.txt")
	return outdir, temp_dir, file_list


def import_file_list(file_list):
	import_tuple = (False,False,False,False,False)
	prev_paired = ""
	paired_separate_list = []
	paired_interleaved_list = []
	unpaired_list = []
	ref_fasta_list = []
	ref_db_list = []
	paired_all_list = []
	with open(file_list) as file_handle:
		for line in file_handle:
			if line.startswith("#") or line == "\n":
				continue
			elif line == "[paired_separate_files]\n":
				print("Importing paired separate files...")
				import_tuple = (True,False,False,False,False)
				continue
			elif line == "[paired_interleaved]\n":
				print("Importing paired interleaved files...")
				import_tuple = (False,True,False,False,False)
				continue
			elif line == "[unpaired]\n":
				print("Importing unpaired files...")
				import_tuple = (False,False,True,False,False)
				continue
			elif line == "[reference_fasta]\n"
				print("Importing reference fasta...")
				import_tuple = (False,False,False,True,False)
				continue
			elif line == "[snpEff_reference_database]\n"
				print("Importing snpEff reference database...")
				import_tuple = (False,False,False,False,True)
				continue
			if import_tuple == (True,False,False,False,False):
				if prev_paired == "":
					print("\n")
					prev_paired = line.strip()
				else:
					paired_separate_list.append((prev_paired,line.strip()))
					prev_paired = ""
			# elif import_tuple == (False,True,False,False,False):
				# paired_interleaved_list.append(line.strip())
			# elif import_tuple == (False,False,True,False,False):
				# unpaired_list.append(line.strip())
			# elif import_tuple == (False,False,False,True,False):
				# ref_fasta_list.append(line.strip())
			# elif import_tuple == (False,False,False,False,True):
				# ref_db_list.append(line.strip())
	return paired_all_list, paired_separate_list, paired_interleaved_list, unpaired_list, ref_fasta_list, ref_db_list


def copy_reference(base_path,ref_fasta_list,temp_dir):
	picard_path = os.path.join(base_path,"picard-tools/picard.jar")
	print("Importing reference sequence: %s" % os.path.basename(ref_fasta_list[0])
	ref_path = os.path.join(temp_dir,"reference.fa")
	with cd(temp_dir):
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"CreateSequenceDictionary","REFERENCE=reference.fa","OUTPUT=reference.dict"])
	shutil.copyfile(ref_fasta_list[0],ref_path)


def copy_ref_db_path(ref_db_list,temp_dir):
	with open(os.path.join(temp_dir,"snpEff_ref_db.txt"),'w') as outfile:
		outfile.write(ref_db_list[0])


def output_keylists(base_path, temp_dir, paired_all_list, unpaired_list):
	for reads1, reads2 in paired_separate_list:
		read1_basename = os.path.basename(reads1)
		read1_filename, ext = os.path.splitext(read1_basename)
		sample_name = read1_filename[:-1].rstrip("._-")
		sample_outfile = "%s.txt" % sample_name
		paired_output1 = "%s.1.paired_trimmed.fastq" % sample_name
		paired_output2 = "%s.2.paired_trimmed.fastq" % sample_name
		unpaired_output1 = "%s.1.unpaired_trimmed.fastq" % sample_name
		unpaired_output2 = "%s.2.unpaired_trimmed.fastq" % sample_name
		with open(os.path.join(base_path,"keys",sample_outfile),'w') as outfile:
			outfile.write("paired\nsplit\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n" % (sample_name,reads1,reads2,paired_output1,paired_output2,unpaired_output1,unpaired_output2))
	for reads in paired_interleaved_list:
		sample_name, file_ext = os.path.splitext(os.path.basename(reads))
		sample_outfile = "%s.txt" % sample_name
		paired_output1 = "%s.1.paired_trimmed.fastq" % sample_name
		paired_output2 = "%s.2.paired_trimmed.fastq" % sample_name
		unpaired_output1 = "%s.1.unpaired_trimmed.fastq" % sample_name
		unpaired_output2 = "%s.2.unpaired_trimmed.fastq" % sample_name
		with open(os.path.join(base_path,"keys",sample_outfile),'w') as outfile:
			outfile.write("paired\ninterleaved\n%s\n%s\n%s\n%s\n%s\n%s\n" % (sample_name,reads,paired_output1,paired_output2,unpaired_output1,unpaired_output2))
	for reads in unpaired_list:
		sample_name, ext = os.path.splitext(os.path.basename(reads))
		sample_outfile = "%s.txt" % sample_name
		output = "%s.trimmed.fastq" % sample_name
		with open(os.path.join(base_path,"keys",sample_outfile),'w') as outfile:
			outfile.write("unpaired\n%s\n%s\n%s\n" % (sample_name,reads,output))


outdir, temp_dir, file_list = setup_directories_and_inputs(base_path)
paired_all_list, paired_separate_list, paired_interleaved_list, unpaired_list, ref_fasta_list, ref_db_list = import_file_list(file_list)
paired_all_list = copy_paired_separate(paired_all_list,paired_separate_list,temp_dir)
paired_all_list = output_split_list(base_path, paired_all_list,paired_interleaved_list,temp_dir)
unpaired_list = copy_unpaired(unpaired_list,temp_dir)
null_output = copy_reference(ref_fasta_list,temp_dir)
null_output = copy_ref_db_path(ref_db_list,temp_dir)
output_keylists(base_path, temp_dir, paired_all_list, unpaired_list)

