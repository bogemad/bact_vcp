#!/usr/bin/env python

import os, sys, shutil, subprocess, collections, gzip, re, vcf, numpy, csv, shlex
from Bio import SeqIO
from contextlib import contextmanager
from multiprocessing import Pool
import pdb


class FASTQRead(collections.namedtuple('FASTQRead', 'name comment seq qual')):
	@property
	def fullname(self):
		if self.comment:
			return '%s %s' % (self.name, self.comment)
		else:
			return self.name

	def __repr__(self):
		if self.comment:
			return '@%s %s\n%s\n+\n%s\n' % (self.name, self.comment, self.seq, self.qual)
		else:
			return '@%s\n%s\n+\n%s\n' % (self.name, self.seq, self.qual)

	def subseq(self, start, end, comment=None):
		if self.comment:
			comment = '%s %s' % (self.comment, comment)

		return FASTQRead(self.name, comment, self.seq[start:end], self.qual[start:end])

	def clone(self, name=None, comment=None, seq=None, qual=None):
		n = name if name else self.name
		c = comment if comment else self.comment
		s = seq if seq else self.seq
		q = qual if qual else self.qual

		return FASTQRead(n, c, s, q)

	def write(self, out):
		out.write(repr(self))


def fastq_read_file(fileobj):
	name = fileobj.next().strip()[1:]

	spl = re.split(r'[ \t]', name, maxsplit=1)
	name = spl[0]
	if len(spl) > 1:
		comment = spl[1]
	else:
		comment = ''

	seq = fileobj.next().strip()
	fileobj.next()
	qual = fileobj.next().strip()

	return FASTQRead(name, comment, seq, qual)


class FASTQ(object):
	def __init__(self, fname=None, fileobj=None):
		self.fname = fname
		self._is_paired = None
		self._is_colorspace = None

		if fileobj:
			self.fileobj = fileobj
		elif fname:
			if fname == '-':
				self.fileobj = sys.stdin
			elif fname[-3:] == '.gz' or fname[-4:] == '.bgz':
				self.fileobj = gzip.open(os.path.expanduser(fname))
			else:
				self.fileobj = open(os.path.expanduser(fname))
		else:
			raise ValueError("Must pass either a fileobj or fname!")

	def tell(self):
		# always relative to uncompressed...
		return self.fileobj.tell()

	def seek(self, pos, whence=0):
		self.fileobj.seek(pos, whence)

	def fetch(self, quiet=False, callback=None):
		# if self.fname and not quiet:
			# eta = ETA(os.stat(self.fname).st_size, fileobj=self.fileobj)
		# else:
			# eta = None

		while True:
			try:
				read = fastq_read_file(self.fileobj)
				# if eta:
					# if callback:
						# eta.print_status(extra=callback())
					# else:
						# eta.print_status(extra=read.name)
				yield read

			except StopIteration:
				break

		# if eta:
			# eta.done()

	def close(self):
		if self.fileobj != sys.stdout:
			self.fileobj.close()

	def check_qualtype(self, num_to_check=10000):
		'''
		Checks a FASTQ file's quality score to see what encoding/scaling is used:
		Sanger, Solexa, or Illumina

		returns "Sanger", "Solexa", "Illumina", or "Unknown"
		'''

		pos = self.tell()

		# these are the differential values, unscaled from chr()
		sanger = (33, 74)  # default sanger is 0->40, but some newer illumina on this scale is 0->41
		solexa = (59, 104)
		illumina = (64, 104)

		sanger_count = 0
		solexa_count = 0
		illumina_count = 0
		unknown_count = 0

		checked = 0
		for read in self.fetch(quiet=True):
			if checked > num_to_check:
				break
			qmax = None
			qmin = None
			for q in [ord(x) for x in read.qual]:
				if qmin is None or q < qmin:
					qmin = q
				if qmax is None or q > qmax:
					qmax = q

			if sanger[0] <= qmin <= qmax <= sanger[1]:
				sanger_count += 1
			elif illumina[0] <= qmin <= qmax <= illumina[1]:
				illumina_count += 1
			elif solexa[0] <= qmin <= qmax <= solexa[1]:
				solexa_count += 1
			else:
				unknown_count += 1
			checked += 1

		self.seek(pos)

		if unknown_count > 0:
			return 'Unknown'  # We don't have any idea about at least one of these reads

		if solexa_count > 0:
			# If there are any reads that fall in the Solexa range,
			# this must be a Solexa scale file. This should be rare.
			return 'Solexa'

		if sanger_count > illumina_count:
			return 'Sanger'
		return 'Illumina'

	@property
	def is_colorspace(self):
		'''
		This works by scanning the first 10 reads that have sequences (aren't Ns
		or 4s). If there are any colorspace values, the entire file is called as
		colorspace.

		It's a bit overkill...
		'''

		if self._is_colorspace is not None:
			return self._is_colorspace

		pos = self.tell()
		self.seek(0)
		self._is_colorspace = None

		valid_basespace = "atcgATCG"
		valid_colorspace = "0123456"

		for read in self.fetch(quiet=True):
			if len(read.seq) < 2:
				continue

			for base in read.seq[1:]:  # skip the first base, in case there is a linker prefix
				if base in valid_colorspace:
					self._is_colorspace = True
					break
				elif base in valid_basespace:
					self._is_colorspace = False
					break
			if self._is_colorspace is not None:
				break

		self.seek(pos)
		return self._is_colorspace

	@property
	def pair_count(self):
		if self.is_paired:  # this actually does the calculation and is cached
			return self._is_paired
		return 1

	@property
	def is_paired(self):
		'''
		Determines if a FASTQ file has paired reads. This returns True is the file has
		paired reads with the same name in consecutive order.
		'''

		if self._is_paired is not None:
			return (self._is_paired > 1)

		pos = self.tell()
		self.seek(0)
		last_name = None
		count = 0

		for read in self.fetch(quiet=True):
			name = read.name.split()[0]
			if last_name:
				if name == last_name:
					count += 1
				else:
					self._is_paired = count
					self.seek(pos)

					return (self._is_paired > 1)
			else:
				last_name = name
				count = 1

		# if there are only 2 reads...
		self._is_paired = count
		self.seek(pos)
		return (self._is_paired > 1)


def unmerge(combined_fname, out_template, gz=False):
	outs = []
	if gz:
		outs.append(gzip.open('%s.1.fastq.gz' % out_template, 'w'))
	else:
		outs.append(open('%s.1.fastq' % out_template, 'w'))

	outidx = 1

	last_read = None
	fq = FASTQ(combined_fname)
	for read in fq.fetch():
		if last_read and last_read.name == read.name:
			outidx += 1
			if len(outs) < outidx:
				if gz:
					outs.append(gzip.open('%s.%s.fastq.gz' % (out_template, outidx), 'w'))
				else:
					outs.append(open('%s.%s.fastq' % (out_template, outidx), 'w'))
			read.write(outs[outidx - 1])
		else:
			outidx = 1
			read.write(outs[0])

		last_read = read

	fq.close()
	for out in outs:
		out.close()


def generator_fetch(generators):
	while True:
		try:
			yield [gen.next() for gen in generators]
		except StopIteration:
			return


def merge(fastqs, split_slashes=False, out=sys.stdout, quiet=False):
	for reads in generator_fetch([fq.fetch(quiet=quiet if i == 0 else False) for i, fq in enumerate(fastqs)]):
		cur_name = None
		for read in reads:
			name = read.name
			comment = read.comment

			if split_slashes and '/' in name:
				spl = name.split('/', 1)
				name = spl[0]
				if read.comment:
					comment = '/%s %s' % (spl[1], read.comment)
				else:
					comment = '/%s' % spl[1]

			if not cur_name:
				cur_name = name
			else:
				if name != cur_name:
					raise ValueError('Files are not paired! Expected: "%s", got "%s"!' % (cur_name, name))

			read.clone(name=name, comment=comment).write(out)


@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)


def replace_dir(dir):
	if os.path.exists(dir):
		shutil.rmtree(dir)
	os.mkdir(dir)


def run_command(command, outputfile):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	while True:
		output = process.stdout.readline()
		if output == '' and process.poll() is not None:
			break
		if output:
			with open(outputfile,'a+') as outhandle:
				outhandle.write(output)
	if process.poll() == 1:
		print "\n\n###########    Error    ###########\n\nThe following command failed to complete successfully:\n\n%s\n\nOutput of this error can be found in:\n\n%s\n\n" % (" ".join(command), outputfile)
		sys.exit(1)


def run_po_command(command, outputfile, errorfile):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = stdout, stderr = process.communicate()
	if output:
		with open(outputfile,'w') as outhandle:
			outhandle.write(output)
	if error:
		with open(errorfile,'a+') as errorhandle:
			errorhandle.write(error)
	if process.returncode == 1:
		print "\n\n###########    Error    ###########\n\nThe following command failed to complete successfully:\n\n%s\n\nOutput of this error can be found in:\n\n%s\n\n" % (" ".join(command), errorfile)
		sys.exit(1)


def setup_directories_and_inputs(base_path, outdir):
	temp_dir = os.path.join(base_path,"intermediate_files")
	replace_dir(temp_dir)
	print "Beginning bact_vcp..."
	if sys.argv[2] == 'true':
		print "Overwriting any previous results."
		replace_dir(outdir)
	elif sys.argv[2] == 'false':
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


def import_file_list(file_list,outdir):
	reads_list = []
	temp_dir = os.path.join(base_path,"intermediate_files")
	with open(file_list) as file_handle:
		for line in file_handle:
			line = re.sub(r'^"#.*"',"",line) # fixes quotes that MS excel sometimes adds to csv files
			if line.startswith("#") or line.startswith(",") or line.strip() == "":
				continue
			data = line.rstrip().split(',')
			reads_list.append((data[0],data[1],data[2],data[3],data[4]))
	return reads_list


def process_reference(base_path,ref_gb,temp_dir,outdir):
	picard_path = os.path.join(base_path, "bin/picard-tools/picard.jar")
	print "Importing reference sequence: %s" % os.path.basename(ref_gb)
	ref_fasta = os.path.join(temp_dir,"reference.fa")
	count = 0
	with open(ref_gb) as original, open(ref_fasta, 'w') as corrected:
		records = SeqIO.parse(original,"gb")
		for record in records:
			record.id = record.name
			record.description = record.name
			count += SeqIO.write(record,corrected,"fasta")
	print "%s sequences written to new fasta file: reference.fa" % count
	outfile = os.path.join(outdir,"Reference_processing_output.txt")
	with cd(temp_dir):
		if os.path.isfile("reference.dict"): os.remove("reference.dict")
		print "Indexing reference sequences."
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"CreateSequenceDictionary","REFERENCE=reference.fa","OUTPUT=reference.dict"], outfile)
		run_command(["smalt","index","-k","13","-s","8","reference","reference.fa"], outfile)
		run_command(["samtools","faidx","reference.fa"], outfile)
	shutil.copy(os.path.join(temp_dir,"reference.fa"),os.path.join(outdir))


def trim_reads(read_data,base_path,temp_dir,outdir,sample_temp_dir,log,results_dir):
	trimo_path = os.path.join(base_path,"bin/trimmomatic/trimmomatic-0.36.jar")
	reads, rgid, sample_name, rgpl, rglb = read_data
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
		run_command(["java",
		"-Xmx2g","-XX:+UseSerialGC",
		"-jar",trimo_path,"PE","-phred33",
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
	shutil.copy(os.path.join(sample_temp_dir,merged_trimmed),os.path.join(results_dir,merged_trimmed))
	shutil.rmtree(sample_temp_dir)
	print "%s: Read trimming complete. Saved to: %s" % (sample_name, os.path.join(results_dir,merged_trimmed))


def sam_to_bam(sam,read_data,base_path,log):
	reads, rgid, sample_name, rgpl, rglb = read_data
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	bam = "%s.raw_alignment.bam" % sample_name
	pre_bam = "%s.raw_no_rg.bam" % sample_name
	run_po_command(["samtools", "view", "-bh", "-q", "10", sam], pre_bam, log)
	run_command(["java",
	"-Xmx2g","-XX:+UseSerialGC",
	"-jar",picard_path,
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
	alignment_quality_check(bam,picard_path,sample_name,log)
	return bam


def alignment_quality_check(bam,picard_path,sample_name,log):
	stats_dir = "alignment_quality_stats"
	os.mkdir(stats_dir)
	prefix = os.path.join(os.path.abspath(stats_dir),sample_name)
	print "%s: Generating alignment quality plots" % sample_name
	run_command(["java",
	"-Xmx2g","-XX:+UseSerialGC",
	"-jar",picard_path,
	"CollectGcBiasMetrics",
	"R=../reference.fa",
	"I=%s" % bam,
	"O=%s_GCBias.txt" % prefix,
	"CHART=%s_GCBias.pdf" % prefix,
	"S=%s_summary_metrics.txt" % prefix,
	"ASSUME_SORTED=true"],log)
	run_command(["java",
	"-Xmx2g","-XX:+UseSerialGC",
	"-jar",picard_path,
	"MeanQualityByCycle",
	"R=../reference.fa",
	"I=%s" % bam,
	"O=%s_Qcycle.txt" % prefix,
	"CHART=%s_Qcycle.pdf" % prefix],log)
	run_command(["java",
	"-Xmx2g","-XX:+UseSerialGC",
	"-jar",picard_path,
	"QualityScoreDistribution",
	"R=../reference.fa",
	"I=%s" % bam,
	"O=%s_Qdist.txt" % prefix,
	"CHART=%s_Qdist.pdf" % prefix],log)
	print "%s: Alignment quality plots complete. Saved to: %s" % (sample_name, os.path.abspath(stats_dir))


def run_smalt(read_data,base_path,temp_dir,outdir,merged_trimmed,sample_temp_dir,log,results_dir):
	#pdb.set_trace()
	reads, rgid, sample_name, rgpl, rglb = read_data
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
	shutil.copytree(os.path.join(sample_temp_dir,"alignment_quality_stats"),os.path.join(results_dir,"alignment_quality_stats"))
	shutil.copyfile(os.path.join(sample_temp_dir,bam),os.path.join(results_dir,bam))
	shutil.rmtree(sample_temp_dir)
	print "%s: Raw alignment bam file complete. Saved to: %s" % (sample_name, os.path.join(results_dir,bam))


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



def process_bam(read_data, base_path,temp_dir,outdir,results_bam,sample_temp_dir,log,results_dir):
	print "%s: Removing PCR duplicates" % read_data[2]
	reads, rgid, sample_name, rgpl, rglb = read_data
	replace_dir(sample_temp_dir)
	bam = os.path.join(sample_temp_dir,"%s.raw_alignment.bam" % sample_name)
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
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
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",markd_input],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"MarkDuplicates",markd_input,markd_output,markd_metrics_arg],log)
		print "%s: PCR duplicates removed" % sample_name
		print "%s: Realigning indels" % sample_name
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","RealignerTargetCreator","-R","../reference.fa","-I",dedup_bam,"-o",realign_intervals],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","IndelRealigner","-R","../reference.fa","-I",dedup_bam,"-targetIntervals",realign_intervals,"-o",preprocessed_reads],log)
		print "%s: Indel realignment complete." % sample_name
		print "%s: Generating depth of coverage analysis" % sample_name
		cov_stats_dir = "coverage_statistics"
		os.mkdir(cov_stats_dir)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","DepthOfCoverage","-R","../reference.fa","-o",os.path.join(cov_stats_dir,cov_stats_name),"-I",preprocessed_reads,"-ct","5","-ct","10"],log)
		print "%s: Depth of coverage analysis complete. Saved to: %s" % (sample_name, os.path.join(results_dir,cov_stats_dir))
		print "%s: Processing of bam alignment complete. Saved to: %s" % (sample_name, os.path.join(results_dir,preprocessed_reads))
	shutil.copytree(os.path.join(sample_temp_dir,cov_stats_dir),os.path.join(results_dir,cov_stats_dir))
	shutil.copy(os.path.join(sample_temp_dir,markd_metrics),os.path.join(results_dir))
	shutil.copy(os.path.join(sample_temp_dir,preprocessed_reads),os.path.join(results_dir))
	shutil.rmtree(sample_temp_dir)



def call_variants(read_data,base_path,temp_dir,outdir,results_preprocessed_reads,sample_temp_dir,log,results_dir):
	print "%s: Calling variants" % read_data[2]
	reads, rgid, sample_name, rgpl, rglb = read_data
	replace_dir(sample_temp_dir)
	preprocessed_reads = "%s.preprocessed_reads.bam" % sample_name
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	raw_variants = "%s.raw_variants.g.vcf" % sample_name
	bbiinput = "INPUT=%s" % preprocessed_reads
	with cd(sample_temp_dir):
		shutil.copyfile(results_preprocessed_reads,preprocessed_reads)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput],log)
		run_command(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","HaplotypeCaller","-R","../reference.fa","-I",preprocessed_reads,"--genotyping_mode","DISCOVERY","-stand_call_conf","30","-stand_emit_conf","10","-ERC", "GVCF", "-o",raw_variants,"--sample_ploidy","1"],log)
	shutil.copyfile(os.path.join(sample_temp_dir,raw_variants),os.path.join(results_dir,raw_variants))
	shutil.rmtree(sample_temp_dir)
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
	snpeff_path = os.path.join(base_path,"bin/snpEff/snpEff.jar")
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
		run_po_command(['java',
		'-Xmx2g', '-XX:+UseSerialGC',
		'-jar', snpeff_path, 
		'-v', 'reference',
		filtered_variants], annotated_variants, log)
	shutil.copyfile(os.path.join(temp_dir,filtered_variants),os.path.join(outdir,filtered_variants))
	shutil.copyfile(os.path.join(temp_dir,annotated_variants),os.path.join(outdir,annotated_variants))
	print "Variant filtering complete. Combined vcf saved to: %s" % os.path.join(outdir,filtered_variants)
	print "Variant annotation complete. Combined vcf saved to: %s" % os.path.join(outdir,annotated_variants)




def run_pipeline(read_data,base_path,temp_dir,outdir):
	results_dir = os.path.join(outdir,read_data[2])
	if os.path.exists(results_dir) == False:
		os.mkdir(results_dir)
	sample_temp_dir = os.path.join(temp_dir,read_data[2])
	log = "%s/%s.log" % (results_dir,read_data[2])
	merged_trimmed = os.path.join(results_dir,"%s.trimmed.fastq.gz" % read_data[2])
	bam = os.path.join(results_dir,"%s.raw_alignment.bam" % read_data[2])
	dedup_bam = os.path.join(results_dir,"%s.dedup_reads.bam" % read_data[2])
	preprocessed_reads = os.path.join(results_dir,"%s.preprocessed_reads.bam" % read_data[2])
	raw_variants = os.path.join(results_dir,"%s.raw_variants.g.vcf" % read_data[2])
	if os.path.exists(merged_trimmed) == False:
		trim_reads(read_data,base_path,temp_dir,outdir,sample_temp_dir,log,results_dir)
	if os.path.exists(bam) == False:
		run_smalt(read_data,base_path,temp_dir,outdir,merged_trimmed,sample_temp_dir,log,results_dir)
	if os.path.exists(preprocessed_reads) == False:
		process_bam(read_data,base_path,temp_dir,outdir,bam,sample_temp_dir,log,results_dir)
	if os.path.exists(raw_variants) == False:
		call_variants(read_data,base_path,temp_dir,outdir,preprocessed_reads,sample_temp_dir,log,results_dir)


def multiprocessing(read_data):
	outdir = os.path.abspath(sys.argv[4])
	base_path = os.path.dirname(os.path.dirname(__file__))
	temp_dir = os.path.join(base_path,"intermediate_files")
	run_pipeline(read_data,base_path,temp_dir,outdir)


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


def skip_line(inhandle,outhandle):
	line = inhandle.next()
	outhandle.write(line)


def add_text_to_line(inhandle,outhandle,text):
	line = inhandle.next().rstrip("\n")
	line += text
	line += "\n"
	outhandle.write(line)




def make_hpc_script(job, base_path, sample_data, job_number, sample_name):
	script_dir = os.path.join(base_path,"scripts")
	job_script_name = "hpc.%s.%s.sh" % (sample_name, job)
	reads, rgid, rgpl, rglb, results_dir, sample_temp_dir, log, merged_trimmed, bam, dedup_bam, preprocessed_reads, raw_variants = sample_data
	hpc_output = "%s.%s.hpc_commandline_output.log" % (job,sample_name)
	hpc_error = "%s.%s.hpc_commandline_error.log" % (job,sample_name)
	hpc_output_dir = os.path.join(results_dir,"commandline_output")
	cmd = '$workdir/scripts/%s.py %s %s %s %s %s %s %s %s %s %s %s %s %s %s' % (job, base_path, sample_name, reads, rgid, rgpl, rglb, results_dir, sample_temp_dir, log, merged_trimmed, bam, dedup_bam, preprocessed_reads, raw_variants)
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
	run_command(["qsub",script])


def del_hpc_script(base_path, sample_name, job):
	script_dir = os.path.join(base_path,"scripts")
	job_script_name = "hpc.%s.%s.sh" % (sample_name, job)
	if os.path.isfile(os.path.join(script_dir,job_script_name)):
		os.remove(os.path.join(script_dir,job_script_name))


def hpc_multiprocessing(base_path, numbered_sample_list, sample_data_dict):
	for job_number, sample_name in numbered_sample_list: #trim_reads loop
		sample_data = sample_data_dict[sample_name]
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
				return


def generate_sample_data_dict(read_list):
	dict = {}
	for read_data in read_list:
		sample_name = read_data[2]
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
	sample_data_dict = generate_sample_data_dict(reads_list)
	numbered_sample_list = [(job_number, sample_name) for (job_number, sample_name) in enumerate(sample_data_dict.keys(),start=1)]
	hpc_multiprocessing(base_path, numbered_sample_list, sample_data_dict)


base_path = os.path.dirname(os.path.dirname(__file__))
file_list = os.path.abspath(sys.argv[1])
ref_gb = os.path.abspath(sys.argv[3])
outdir = os.path.abspath(sys.argv[4])



if __name__ == "__main__":
	temp_dir = setup_directories_and_inputs(base_path, outdir)
	reads_list = import_file_list(file_list,outdir)
	null_output = process_reference(base_path,ref_gb,temp_dir,outdir)
	setup_snpEff(ref_gb,base_path,outdir)
	try:
		if sys.argv[5] == 'true':
			run_pipeline_hpc(reads_list,base_path,temp_dir,outdir)
		if sys.argv[5] == 'false':
			p = Pool(int(sys.argv[6]))
			p.map_async(multiprocessing,reads_list).get(9999999)
			p.close()
			p.join()
		if os.path.exists(os.path.join(outdir,"raw_variants.vcf")) == False:
			genotype_variants(reads_list,base_path,temp_dir,outdir)
		if os.path.exists(os.path.join(outdir,"annotated_variants.vcf")) == False:
			annotate_variants(base_path,temp_dir,outdir)
	except:
		shutil.rmtree(temp_dir)
		reset_snpEff(ref_gb,base_path)
		raise
	else:
		shutil.rmtree(temp_dir)
		reset_snpEff(ref_gb,base_path)
		print "All jobs completed successfully."


