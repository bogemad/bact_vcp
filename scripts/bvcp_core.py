#!/usr/bin/env python

import os, gzip, shutil
from bvcp_util import replace_dir, cd, copy_to_results_dir, run_command, run_po_command 
from fastqutils import unmerge, merge, FASTQ

def trim_reads(read_data,base_path,temp_dir,sample_temp_dir,log,results_dir):
	reads, rgid, sample_name, rgpl, rglb, outdir = read_data
	replace_dir(sample_temp_dir)
	paired_output1 = "%s.trimmed_paired1.fastq" % sample_name
	paired_output2 = "%s.trimmed_paired2.fastq" % sample_name
	unpaired_output1 = "%s.trimmed_unpaired1.fastq" % sample_name
	unpaired_output2 = "%s.trimmed_unpaired2.fastq" % sample_name
	merged_trimmed = "%s.trimmed.fastq.gz" % sample_name
	print "%s: Importing reads" % sample_name 
	trim_log = "%s.trimlog" % sample_name
	reads1 = "%s.1.fastq" % sample_name
	reads2 = "%s.2.fastq" % sample_name
	log = os.path.abspath("%s/%s.log" % (results_dir,sample_name))
	with cd(sample_temp_dir):
		print "%s: Processing reads file for trimming" % sample_name 
		unmerge(reads, sample_name, gz=False)
		print "%s: Trimming reads" % sample_name 
		cut_adapt_out1 = "ca_%s" % reads1
		cut_adapt_out2 = "ca_%s" % reads2
		run_command(["cutadapt",
		"-g","AGATGTGTATAAGAGACAG",
		"-a","CTGTCTCTTATACACATCT",
		"-g","AGATGTGTATAAGAGACAG",
		"-a","CTGTCTCTTATACACATCT",
		"-g","TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
		"-a","CTGTCTCTTATACACATCTGACGCTGCCGACGA",
		"-g","GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
		"-a","CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
		"-o",cut_adapt_out1,"-p",cut_adapt_out2,
		reads1, reads2],log)
		run_command(["trimmomatic",
		"-Xmx2g","-XX:+UseSerialGC",
		"PE","-phred33",
		"-threads","1",
		"-trimlog",trim_log,
		cut_adapt_out1,cut_adapt_out2,
		paired_output1,unpaired_output1,
		paired_output2,unpaired_output2,
		"LEADING:20","TRAILING:20",
		"SLIDINGWINDOW:4:15","MINLEN:70"],log)
		print "%s: Merging trimmed reads" % sample_name
		with gzip.open(merged_trimmed, 'wb') as trimmed_merged:
			fastqs = [FASTQ(x) for x in (paired_output1,paired_output2)]
			merge(fastqs, split_slashes=False, out=trimmed_merged, quiet=False)
	copy_to_results_dir(os.path.join(sample_temp_dir,merged_trimmed), results_dir)
	print "%s: Read trimming complete. Saved to: %s" % (sample_name, os.path.join(results_dir,merged_trimmed))


def alignment_quality_check(bam,sample_name,log):
	stats_dir = "alignment_quality_stats"
	os.mkdir(stats_dir)
	prefix = os.path.join(os.path.abspath(stats_dir),sample_name)
	print "%s: Generating alignment quality plots" % sample_name
	run_command(["picard",
	"-Xmx2g","-XX:+UseSerialGC",
	"CollectGcBiasMetrics",
	"R=../reference.fa",
	"I=%s" % bam,
	"O=%s_GCBias.txt" % prefix,
	"CHART=%s_GCBias.pdf" % prefix,
	"S=%s_summary_metrics.txt" % prefix,
	"ASSUME_SORTED=true"],log)
	run_command(["picard",
	"-Xmx2g","-XX:+UseSerialGC",
	"MeanQualityByCycle",
	"R=../reference.fa",
	"I=%s" % bam,
	"O=%s_Qcycle.txt" % prefix,
	"CHART=%s_Qcycle.pdf" % prefix],log)
	run_command(["picard",
	"-Xmx2g","-XX:+UseSerialGC",
	"QualityScoreDistribution",
	"R=../reference.fa",
	"I=%s" % bam,
	"O=%s_Qdist.txt" % prefix,
	"CHART=%s_Qdist.pdf" % prefix],log)
	print "%s: Alignment quality plots complete. Saved to: %s" % (sample_name, os.path.abspath(stats_dir))



def sam_to_bam(sam,read_data,base_path,log):
	reads, rgid, sample_name, rgpl, rglb, outdir = read_data
	bam = "%s.raw_alignment.bam" % sample_name
	pre_bam = "%s.raw_no_rg.bam" % sample_name
	run_po_command(["samtools", "view", "-bh", "-q", "10", sam], pre_bam, log)
	run_command(["picard",
	"-Xmx2g","-XX:+UseSerialGC",
	"AddOrReplaceReadGroups",
	"I=%s" % pre_bam,
	"O=%s" % bam,
	"RGID=%s" % rgid,
	"RGSM=%s" % sample_name,
	"RGPL=%s" % rgpl.lower(),
	"RGLB=%s" % rglb,
	"RGPU=NA",
	"SORT_ORDER=coordinate"],log)
	os.remove(sam)
	os.remove(pre_bam)
	alignment_quality_check(bam,sample_name,log)
	return bam


def run_smalt(read_data,base_path,temp_dir,merged_trimmed,sample_temp_dir,log,results_dir):
	reads, rgid, sample_name, rgpl, rglb, outdir = read_data
	replace_dir(sample_temp_dir)
	reads1 = "%s.1.fastq" % sample_name
	reads2 = "%s.2.fastq" % sample_name
	with cd(sample_temp_dir):
		print "%s: Processing reads for mapping" % sample_name
		unmerge(merged_trimmed, sample_name, gz=False)
		sam = "%s.sam" % sample_name
		print "%s: Mapping reads" % sample_name
		run_command(["smalt","map","-i","1500","-n","1","-o",sam,"../reference",reads1,reads2],log)
		print "%s: Read mapping complete" % sample_name
		print "%s: Generating raw alignment bam file" % sample_name
		bam = sam_to_bam(sam,read_data,base_path,log)
	copy_to_results_dir(os.path.join(sample_temp_dir,"alignment_quality_stats"), results_dir)
	copy_to_results_dir(os.path.join(sample_temp_dir,bam), results_dir)
	print "%s: Raw alignment bam file complete. Saved to: %s" % (sample_name, os.path.join(results_dir,bam))


def process_bam(read_data, base_path,temp_dir,results_bam,sample_temp_dir,log,results_dir):
	print "%s: Removing PCR duplicates" % read_data[2]
	reads, rgid, sample_name, rgpl, rglb, outdir = read_data
	replace_dir(sample_temp_dir)
	bam = os.path.join(sample_temp_dir,"%s.raw_alignment.bam" % sample_name)
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	dedup_bam = "%s.dedup_reads.bam" % sample_name
	markd_input = "INPUT=%s" % bam
	markd_output = "OUTPUT=%s" % dedup_bam
	markd_metrics = "%s.picard_md_metrics.txt" % sample_name
	markd_metrics_arg = "METRICS_FILE=%s" % markd_metrics
	realign_intervals = "%s.target.intervals" % sample_name
	preprocessed_reads = "%s.preprocessed_reads.bam" % sample_name
	bbiinput = "INPUT=%s" % dedup_bam
	cov_stats_name = "%s.coverage_stats" % sample_name
	with cd(sample_temp_dir):
		shutil.copyfile(results_bam,bam)
		run_command(["picard","-Xmx2g","-XX:+UseSerialGC","BuildBamIndex",markd_input],log)
		run_command(["picard","-Xmx2g","-XX:+UseSerialGC","MarkDuplicates",markd_input,markd_output,markd_metrics_arg],log)
		print "%s: PCR duplicates removed" % sample_name
		print "%s: Realigning indels" % sample_name
		run_command(["picard","-Xmx2g","-XX:+UseSerialGC","BuildBamIndex",bbiinput],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","RealignerTargetCreator","-R","../reference.fa","-I",dedup_bam,"-o",realign_intervals],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","IndelRealigner","-R","../reference.fa","-I",dedup_bam,"-targetIntervals",realign_intervals,"-o",preprocessed_reads],log)
		print "%s: Indel realignment complete." % sample_name
		print "%s: Generating depth of coverage analysis" % sample_name
		cov_stats_dir = "coverage_statistics"
		os.mkdir(cov_stats_dir)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","DepthOfCoverage","-R","../reference.fa","-o",os.path.join(cov_stats_dir,cov_stats_name),"-I",preprocessed_reads,"-ct","5","-ct","10"],log)
		print "%s: Depth of coverage analysis complete. Saved to: %s" % (sample_name, os.path.join(results_dir,cov_stats_dir))
		print "%s: Processing of bam alignment complete. Saved to: %s" % (sample_name, os.path.join(results_dir,preprocessed_reads))
	copy_to_results_dir(os.path.join(sample_temp_dir,cov_stats_dir), results_dir)
	copy_to_results_dir(os.path.join(sample_temp_dir,markd_metrics), results_dir)
	copy_to_results_dir(os.path.join(sample_temp_dir,preprocessed_reads), results_dir)


def call_variants(read_data,base_path,temp_dir,results_preprocessed_reads,sample_temp_dir,log,results_dir):
	print "%s: Calling variants" % read_data[2]
	reads, rgid, sample_name, rgpl, rglb, outdir = read_data
	replace_dir(sample_temp_dir)
	preprocessed_reads = "%s.preprocessed_reads.bam" % sample_name
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	raw_variants = "%s.raw_variants.g.vcf" % sample_name
	bbiinput = "INPUT=%s" % preprocessed_reads
	with cd(sample_temp_dir):
		shutil.copyfile(results_preprocessed_reads,preprocessed_reads)
		run_command(["picard","-Xmx2g","-XX:+UseSerialGC","BuildBamIndex",bbiinput],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","HaplotypeCaller","-R","../reference.fa","-I",preprocessed_reads,"--genotyping_mode","DISCOVERY","-stand_call_conf","30","-stand_emit_conf","10","-ERC", "GVCF", "-o",raw_variants,"--sample_ploidy","1"],log)
	copy_to_results_dir(os.path.join(sample_temp_dir,raw_variants), results_dir)
	print "%s: Variant calling complete. Raw variants vcf saved to: %s" % (sample_name, os.path.join(results_dir,raw_variants))


def annotate_variants(base_path,temp_dir,outdir):
	print "First step filtering of variants"
	raw_variants = "raw_variants.vcf"
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	raw_snps = "raw_snps.vcf"
	raw_indels = "raw_indels.vcf"
	filtered_snps = "filtered_snps.vcf"
	filtered_indels = "filtered_indels.vcf"
	filtered_variants = "filtered_variants.vcf"
	annotated_variants = "annotated_variants.vcf"
	log = "%s/variant_filtering_annotation.log" % (outdir)
	with cd(temp_dir):
		shutil.copyfile(os.path.join(outdir,raw_variants),raw_variants)
		run_command(["java",
		"-Xmx2g","-XX:+UseSerialGC",
		"-jar",gatk_path,
		"-T","SelectVariants",
		"-R","reference.fa",
		"-V",raw_variants,
		"-selectType","SNP",
		"-o",raw_snps],log)
		run_command(["java",
		"-Xmx2g","-XX:+UseSerialGC",
		"-jar",gatk_path,
		"-T","SelectVariants",
		"-R","reference.fa",
		"-V",raw_variants,
		"-selectType","INDEL",
		"-o",raw_indels],log)
		run_command(['java',
		'-Xmx2g', '-XX:+UseSerialGC',
		'-jar', gatk_path,
		'-T', 'VariantFiltration',
		'-R', 'reference.fa',
		'-V', raw_snps,
		'--filterExpression', '"QD < 2.0"', '--filterName', '"low_qual_by_depth"',
		'--filterExpression', '"MQ < 40.0"', '--filterName', '"low_RMS_mapping_quality"',
		'--filterExpression', '"FS > 60.0"', '--filterName', '"strand_bias"',
		'--filterExpression', '"SOR > 3.0"', '--filterName', '"high_strand_odds_ratio"',
		'--filterExpression', '"MQRankSum < -12.5"', '--filterName', '"low_MQRankSum"',
		'--filterExpression', '"ReadPosRankSum < -8.0"', '--filterName', '"low_ReadPosRankSum"',
		'-o', filtered_snps], log)
		run_command(['java', 
		'-Xmx2g', '-XX:+UseSerialGC', 
		'-jar', gatk_path,
		'-T', 'VariantFiltration',
		'-R', 'reference.fa',
		'-V', raw_indels,
		'--filterExpression', '"QD < 2.0"', '--filterName', '"low_qual_by_depth"',
		'--filterExpression', '"ReadPosRankSum < -20.0"', '--filterName', '"low_ReadPosRankSum"',
		'--filterExpression', '"FS > 200.0"', '--filterName', '"strand_bias"',
		'--filterExpression', '"SOR > 10.0"', '--filterName', '"high_strand_odds_ratio"',
		'-o', filtered_indels],log)
		print "Annotating variants with snpEff"
		run_command(["java", 
		"-Xmx2g", "-XX:+UseSerialGC",
		"-jar", gatk_path,
		"-T", "CombineVariants",
		"-R", "reference.fa",
		"--variant:snp", filtered_snps,
		"--variant:indel", filtered_indels,
		"-o", filtered_variants,
		"-genotypeMergeOptions", "PRIORITIZE",
		"-priority", "snp,indel"], log)
		run_po_command(['snpEff',
		'-Xmx2g', '-XX:+UseSerialGC',
		'-v', 'reference',
		filtered_variants], annotated_variants, log)
	copy_to_results_dir(os.path.join(temp_dir,filtered_variants), outdir)
	copy_to_results_dir(os.path.join(temp_dir,annotated_variants), outdir)
	print "Variant filtering complete. Combined vcf saved to: %s" % os.path.join(outdir,filtered_variants)
	print "Variant annotation complete. Combined vcf saved to: %s" % os.path.join(outdir,annotated_variants)


def genotype_variants(reads_list,base_path,temp_dir,outdir):
	print "Merging raw variant files"
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	log = "%s/variant_filtering_annotation.log" % (outdir)
	variant_files = []
	for reads_data in reads_list:
		raw_variants = "%s.raw_variants.g.vcf" % reads_data[2]
		variant_files.append(os.path.join(outdir,reads_data[2],raw_variants))
	cmd = ["java", "-Xmx2g", "-XX:+UseSerialGC", "-jar", gatk_path, "-T", "GenotypeGVCFs", "-R", "reference.fa"]
	for raw_variants in variant_files:
		cmd.append("--variant")
		cmd.append(raw_variants)
	cmd.append("-o")
	cmd.append("raw_variants.vcf")
	with cd(temp_dir):
		run_command(cmd, log)
	shutil.copy(os.path.join(temp_dir,"raw_variants.vcf"),outdir)


#Counts number of snps in a 200 bp window across whole genome
# def count_variants_and_plot(filtered_variants):
	# count_list = []
	# win_out = []
	# window = 1000
	# records = [record for record in vcf.Reader(open(filtered_variants, 'r'))]
	# last_entry = records[-1].POS
	# snp_positions = numpy.zeros(last_entry+1, dtype=bool)
	# for record in records:
		# snp_positions[record.POS] = True
	# temp_counts = []
	# for n in xrange(0, last_entry, window):
		# count = numpy.count_nonzero(snp_positions[n:n+window])
		# temp_counts.append(count)
		# if i == 0:
			# win_out.append(n)
	# if i == 0:
		# count_list.append(win_out)	
	# count_list.append(temp_counts)

	# with open("snp_density.csv", "wb") as fh_out:
		# writer = csv.writer(fh_out)
		# files.insert(0, 'position')
		# writer.writerow(files)
		# writer.writerows(zip(*count_list))