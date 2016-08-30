#!/usr/bin/env python

import os, shutil
from bvcp_util import skip_line, add_text_to_line, run_command

def make_hpc_script(job, base_path, sample_data, job_number, sample_name):
	script_dir = os.path.join(base_path,"scripts")
	job_script_name = "hpc.%s.%s.sh" % (sample_name, job)
	reads, rgid, rgpl, rglb, results_dir, sample_temp_dir, log, merged_trimmed, bam, dedup_bam, preprocessed_reads, raw_variants = sample_data
	hpc_output = "%s.%s.hpc_commandline_output.log" % (job,sample_name)
	hpc_error = "%s.%s.hpc_commandline_error.log" % (job,sample_name)
	hpc_output_dir = os.path.join(results_dir,"commandline_output")
	cmd = '$bact_vcp_path/scripts/%s.py %s %s %s %s %s %s %s %s %s %s %s %s %s %s' % (job, base_path, sample_name, reads, rgid, rgpl, rglb, results_dir, sample_temp_dir, log, merged_trimmed, bam, dedup_bam, preprocessed_reads, raw_variants)
	if os.path.isdir(hpc_output_dir) == False:
		os.mkdir(hpc_output_dir)
	with open(os.path.join(script_dir,"hpc.job_script.sh")) as script_template:
		with open(os.path.join(script_dir,job_script_name),'w') as job_script:
			skip_line(script_template,job_script)
			add_text_to_line(script_template,job_script,"%s.%d" % (job, job_number))
			skip_line(script_template,job_script)
			skip_line(script_template,job_script)
			skip_line(script_template,job_script)
			skip_line(script_template,job_script)
			add_text_to_line(script_template,job_script,os.path.join(hpc_output_dir,hpc_output))
			add_text_to_line(script_template,job_script,os.path.join(hpc_output_dir,hpc_error))
			skip_line(script_template,job_script)
			add_text_to_line(script_template,job_script,sample_name)
			add_text_to_line(script_template,job_script,job)
			add_text_to_line(script_template,job_script,base_path)
			skip_line(script_template,job_script)
			skip_line(script_template,job_script)
			skip_line(script_template,job_script)
			skip_line(script_template,job_script)
			add_text_to_line(script_template,job_script,cmd)


def run_hpc_script(base_path, sample_name, job):
	job_script_name = "hpc.%s.%s.sh" % (sample_name, job)
	script = os.path.join(base_path,"scripts",job_script_name)
	run_command(["qsub",script],"/dev/null")


def del_hpc_script(base_path, sample_name, job):
	script_dir = os.path.join(base_path,"scripts")
	job_script_name = "hpc.%s.%s.sh" % (sample_name, job)
	if os.path.isfile(os.path.join(script_dir,job_script_name)):
		os.remove(os.path.join(script_dir,job_script_name))


def multiprocessing_cleanup(sample_data_dict):
	for sample_name in sample_data_dict.keys():
		sample_data = sample_data_dict[sample_name]
		shutil.rmtree(sample_data[5])
		run_command(['conda', 'remove', '--name', sample_name, '--all'], sample_data[6])
		run_command(['conda', 'remove', '--name', sample_name, '--all'], sample_data[6])


def hpc_multiprocessing(base_path, numbered_sample_list, sample_data_dict):
	for job_number, sample_name in numbered_sample_list: #trim_reads loop
		sample_data = sample_data_dict[sample_name]
		run_command(['conda','create','-n',sample_name,'--clone','venv'],sample_data[6])
		if os.path.exists(sample_data[7]) == False: # if trimmed reads don't exist, start trim_reads
			job = "trim_reads"
			make_hpc_script(job, base_path, sample_data, job_number, sample_name)
			run_hpc_script(base_path, sample_name, job)
	smalt_submitted_list = []
	bam_submitted_list = []
	var_submitted_list = []
	job_completed_list = []
	while True:
		for job_number, sample_name in numbered_sample_list: #run_smalt loop
			sample_data = sample_data_dict[sample_name]
			if os.path.exists(sample_data[7]) == False: # if trimmed reads don't exist, don't do anything yet
				continue
			if sample_name not in smalt_submitted_list: #if script hasn't been submitted
				if os.path.exists(sample_data[8]) == False: # if raw bam doesn't exist, start run_smalt
					del_hpc_script(base_path, sample_name, "trim_reads")
					job = "run_smalt"
					smalt_submitted_list.append(sample_name)
					make_hpc_script(job, base_path, sample_data, job_number, sample_name)
					run_hpc_script(base_path, sample_name, job)
				elif os.path.exists(sample_data[8]) == True: # if raw bam already exists, add to job submitted list
					smalt_submitted_list.append(sample_name)
		for job_number, sample_name in numbered_sample_list: #pre_bam loop
			sample_data = sample_data_dict[sample_name]
			if os.path.exists(sample_data[8]) == False:
				continue
			if sample_name not in bam_submitted_list: #if script hasn't been submitted
				if os.path.exists(sample_data[10]) == False:
					del_hpc_script(base_path, sample_name, "run_smalt")
					job = "pre_bam"
					bam_submitted_list.append(sample_name)
					make_hpc_script(job, base_path, sample_data, job_number, sample_name)
					run_hpc_script(base_path, sample_name, job)
				elif os.path.exists(sample_data[10]) == True:
					bam_submitted_list.append(sample_name)
		for job_number, sample_name in numbered_sample_list: #call_var loop
			sample_data = sample_data_dict[sample_name]
			if os.path.exists(sample_data[10]) == False:
				continue
			if sample_name not in var_submitted_list: #if script hasn't been submitted
				if os.path.exists(sample_data[11]) == False:
					del_hpc_script(base_path, sample_name, "pre_bam")
					job = "call_var"
					var_submitted_list.append(sample_name)
					make_hpc_script(job, base_path, sample_data, job_number, sample_name)
					run_hpc_script(base_path, sample_name, job)
				elif os.path.exists(sample_data[11]) == True:
					var_submitted_list.append(sample_name)
		for job_number, sample_name in numbered_sample_list: #wait for all call_var to finish loop
			sample_data = sample_data_dict[sample_name]
			if sample_name not in job_completed_list:
				if os.path.exists(sample_data[11]) == True:
					del_hpc_script(base_path, sample_name, "call_var")
					job_completed_list.append(sample_name)
			if len(numbered_sample_list) == len(job_completed_list):
				multiprocessing_cleanup(sample_data_dict)
				return


def generate_sample_data_dict(read_list, temp_dir):
	dict = {}
	for read_data in read_list:
		sample_name = read_data[2]
		outdir = read_data[5]
		results_dir = os.path.join(outdir,sample_name)
		if os.path.exists(results_dir) == False:
			os.mkdir(results_dir)
		sample_temp_dir = os.path.join(temp_dir,sample_name)
		log = "%s/%s.log" % (results_dir,sample_name)
		merged_trimmed = os.path.join(results_dir,"%s.trimmed.fastq.gz" % sample_name)
		bam = os.path.join(results_dir,"%s.raw_alignment.bam" % sample_name)
		dedup_bam = os.path.join(results_dir,"%s.dedup_reads.bam" % sample_name)
		preprocessed_reads = os.path.join(results_dir,"%s.preprocessed_reads.bam" % sample_name)
		raw_variants = os.path.join(results_dir,"%s.raw_variants.g.vcf" % sample_name)
		dict[sample_name] = read_data[0], read_data[1], read_data[3], read_data[4], results_dir, sample_temp_dir, log, merged_trimmed, bam, dedup_bam, preprocessed_reads, raw_variants
	return dict


def run_pipeline_hpc(reads_list,base_path,temp_dir,outdir):
	sample_data_dict = generate_sample_data_dict(reads_list, temp_dir)
	numbered_sample_list = [(job_number, sample_name) for (job_number, sample_name) in enumerate(sample_data_dict.keys(),start=1)]
	hpc_multiprocessing(base_path, numbered_sample_list, sample_data_dict)


