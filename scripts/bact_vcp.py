#!/usr/bin/env python

import os, sys, shutil
from bvcp_ref import setup_directories_and_inputs, import_file_list, process_reference, setup_snpEff, reset_snpEff
from bvcp_hpc import run_pipeline_hpc
from bvcp_local import run_pipeline_local
from bvcp_core import genotype_variants, annotate_variants


if __name__ == "__main__":
	base_path = os.path.dirname(os.path.dirname(__file__))
	file_list = os.path.abspath(sys.argv[1])
	force = sys.argv[2]
	ref_gb = os.path.abspath(sys.argv[3])
	outdir = os.path.abspath(sys.argv[4])
	temp_dir = setup_directories_and_inputs(base_path, outdir, force)
	reads_list = import_file_list(file_list,outdir,temp_dir)
	null_output = process_reference(base_path,ref_gb,temp_dir,outdir)
	setup_snpEff(ref_gb,base_path,outdir)
	#try:
	if sys.argv[5] == 'true':
		run_pipeline_hpc(reads_list,base_path,temp_dir,outdir)
	elif sys.argv[5] == 'false':
		run_pipeline_local(sys.argv[6], reads_list)
	if os.path.exists(os.path.join(outdir,"raw_variants.vcf")) == False:
		genotype_variants(reads_list,base_path,temp_dir,outdir)
	if os.path.exists(os.path.join(outdir,"annotated_variants.vcf")) == False:
		annotate_variants(base_path,temp_dir,outdir)
	# except:
		# shutil.rmtree(temp_dir)
		# reset_snpEff(ref_gb,base_path)
		# raise
	# else:
		# shutil.rmtree(temp_dir)
		# reset_snpEff(ref_gb,base_path)
		# print "All jobs completed successfully."


