#!/usr/bin/env python

import os, sys, shutil, re
from bvcp_util import replace_dir, run_command, cd
from Bio import SeqIO

def setup_directories_and_inputs(base_path, outdir, force):
	temp_dir = os.path.join(base_path,"intermediate_files")
	replace_dir(temp_dir)
	print "Beginning bact_vcp..."
	if force == 'true':
		print "Overwriting any previous results."
		replace_dir(outdir)
	elif force == 'false':
		if os.path.isdir(outdir) == False:
			print "No previous results found. Creating new results directory..."
			os.mkdir(outdir)
		else:
			print "Previous results found. Adding to current analysis..."
	return temp_dir


def codon_table_lookup(ref_gb):
	codon_table = 0
	codon_table_d = {}
	export_codon_table_d = {}
	codon_table_d[1] = "Standard"
	codon_table_d[2] = "Vertebrate_Mitochondrial"
	codon_table_d[3] = "Yeast_Mitochondrial"
	codon_table_d[4] = "Mycoplasma"
	codon_table_d[5] = "Invertebrate_Mitochondrial"
	codon_table_d[6] = "Ciliate_Nuclear"
	codon_table_d[9] = "Echinoderm_Mitochondrial"
	codon_table_d[10] = "Euplotid_Nuclear"
	codon_table_d[11] = "Bacterial_and_Plant_Plastid"
	codon_table_d[12] = "Alternative_Yeast_Nuclear"
	codon_table_d[13] = "Ascidian_Mitochondrial"
	codon_table_d[14] = "Alternative_Flatworm_Mitochondrial"
	codon_table_d[16] = "Chlorophycean_Mitochondrial"
	codon_table_d[21] = "Trematode_Mitochondrial"
	codon_table_d[22] = "Scenedesmus_obliquus_Mitochondrial"
	codon_table_d[23] = "Thraustochytrium_Mitochondrial"
	try:
		for rec in SeqIO.parse(ref_gb,"gb"):
			for feature in rec.features:
				if feature.type == 'CDS':
					t_l = feature.qualifiers['transl_table']
					if codon_table == 0:
						codon_table = int(t_l[0])
					elif codon_table != int(t_l[0]):
						print "Error: Multiple /transl_table values found in Genbank file, please check file."
						sys.exit(1)
			print "Using %s codon table for chromosome/contig %s" % (codon_table_d[codon_table],rec.name)
			export_codon_table_d[rec.name] = codon_table_d[codon_table]
	except:
		print "Error: /transl_table not found for Genbank file."
		raise
		sys.exit(1)
	return export_codon_table_d


def setup_snpEff(ref_gb,base_path,outdir):
	outfile = os.path.join(outdir,"Reference_processing_output.txt")
	snpeff_path = os.path.join(base_path,"bin/snpEff/snpEff.jar")
	snpeff_dir = os.path.join(base_path,"bin/snpEff")
	snpEff_config = os.path.join(snpeff_dir,"snpEff.config")
	snpEff_config_backup = os.path.join(snpeff_dir,"snpEff.config.bak")
	shutil.copyfile(snpEff_config,snpEff_config_backup)
	data_dir = os.path.join(snpeff_dir,"data")
	reference_data_dir = os.path.join(data_dir,"reference")
	if os.path.isdir(data_dir) == True:
		shutil.rmtree(data_dir)
	os.mkdir(data_dir)
	os.mkdir(reference_data_dir)
	shutil.copyfile(ref_gb,os.path.join(reference_data_dir,"genes.gbk"))
	codon_table_names = codon_table_lookup(ref_gb)
	with open(snpEff_config,'a') as snp_Eff_config_handle:
		snp_Eff_config_handle.write("# reference\n")
		snp_Eff_config_handle.write("reference.genome : reference\n")
		count = 1
		for rec in SeqIO.parse(ref_gb,"gb"):
			if count == 1:
				outline = "reference.chromosomes : %s" % rec.name
			else:
				outline = ", %s" % rec.name
			snp_Eff_config_handle.write(outline)
			count += 1
		snp_Eff_config_handle.write("\n")
		for rec in SeqIO.parse(ref_gb,"gb"):
			outline = "reference.%s.codonTable : %s\n" % (rec.name,codon_table_names[rec.name])
			snp_Eff_config_handle.write(outline)
	run_command(["java","-jar",snpeff_path,"build","-genbank","-v","reference"],outfile)


def reset_snpEff(ref_gb,base_path):
	snpeff_dir = os.path.join(base_path,"bin/snpEff")
	snpEff_config = os.path.join(snpeff_dir,"snpEff.config")
	snpEff_config_backup = os.path.join(snpeff_dir,"snpEff.config.bak")
	data_dir = os.path.join(snpeff_dir,"data")
	os.remove(snpEff_config)
	os.rename(snpEff_config_backup,snpEff_config)
	shutil.rmtree(data_dir)


def import_file_list(file_list,outdir,temp_dir):
	reads_list = []
	with open(file_list) as file_handle:
		for line in file_handle:
			line = re.sub(r'^"#.*"',"",line) # fixes quotes that MS excel sometimes adds to csv files
			if line.startswith("#") or line.startswith(",") or line.strip() == "":
				continue
			data = line.rstrip().split(',')
			reads_list.append((data[0],data[1],data[2],data[3],data[4],outdir))
	return reads_list


def process_reference(base_path,ref_gb,temp_dir,outdir):
	picard_path = os.path.join(base_path, "bin/picard-tools/picard.jar")
	print "Importing reference sequence: %s" % os.path.basename(ref_gb)
	ref_fasta = os.path.join(temp_dir,"reference.fa")
	new_ref_gb = os.path.join(temp_dir,"reference.gbk")
	count = 0
	gbk_count = 0
	with open(ref_gb) as original, open(ref_fasta, 'w') as fasta, open(new_ref_gb, 'w') as gbk:
		records = SeqIO.parse(original,"gb")
		for record in records:
			record.id = record.name
			record.description = record.name
			count += SeqIO.write(record,fasta,"fasta")
			gbk_count += SeqIO.write(record,gbk,"gb")
	print "%s sequences written to new fasta file: reference.fa" % count
	print "%s sequences written to new gbk file: reference.gbk" % gbk_count
	outfile = os.path.join(outdir,"Reference_processing_output.txt")
	with cd(temp_dir):
		if os.path.isfile("reference.dict"): os.remove("reference.dict")
		print "Indexing reference sequences."
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"CreateSequenceDictionary","REFERENCE=reference.fa","OUTPUT=reference.dict"], outfile)
		run_command(["smalt","index","-k","13","-s","8","reference","reference.fa"], outfile)
		run_command(["samtools","faidx","reference.fa"], outfile)
	shutil.copy(os.path.join(temp_dir,"reference.fa"),os.path.join(outdir))

