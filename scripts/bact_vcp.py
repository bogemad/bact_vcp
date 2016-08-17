#!/usr/bin/env python

import os, sys, shutil, subprocess, collections, gzip, re, vcf, numpy, csv
from Bio import SeqIO
from eta import ETA
from contextlib import contextmanager
from multiprocessing import Pool



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


def run_process(cmd):
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	output, err = p.communicate()
	if output:
		print output
	if err:
		print err

def run_shell_process(cmd):
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	output, err = p.communicate()
	if output:
		print output
	if err:
		print err


def run_piped_shell_process(cmd, outputfile):
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	output, err = p.communicate()
	with open(outputfile,'w') as output_handle:
		output_handle.write(output)
	if err:
		print err


def run_piped_shell_process_outerr(cmd, outputfile):
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	output, err = p.communicate()
	with open(outputfile,'w') as output_handle:
		output_handle.write(output)
	if err:
		print err



# def checkpoint_passed(sample_name,checkpoint):
	# temp_dir = os.path.join(base_path,"intermediate_files")
	# check_path = os.path.join(temp_dir,".%s.checkpoint" % sample_name)
	# with open(check_path,'w') as checkpoint_file:
		# checkpoint_file.write(checkpoint)
	# print "%s: %s checkpoint passed!" % (sample_name,checkpoint)


# def checkpoint_check(sample_name):
	# temp_dir = os.path.join(base_path,"intermediate_files")
	# check_path = os.path.join(temp_dir,".%s.checkpoint" % sample_name)
	# with open(check_path) as checkpoint_file:
		# checkpoint = checkpoint_file.next()
	# print "%s: Starting from %s checkpoint" % (sample_name,checkpoint)
	# return checkpoint


def setup_directories_and_inputs(base_path, outdir):
	temp_dir = os.path.join(base_path,"intermediate_files")
	checkpoint_temp = os.path.join(base_path,"checkpoint_files")
	outdir = os.path.abspath(outdir)
	try:
		if sys.argv[2] == 'false':
			if sys.argv[6] == 'true':
				print "Checkpoint run"
				if os.path.isdir(outdir) == False:
					os.mkdir(outdir)
				if os.path.isdir(checkpoint_temp) == True:
					print "Importing checkpoint directory"
					os.rename(checkpoint_temp,temp_dir)
				else:
					if os.path.isdir(temp_dir) == True:
						shutil.rmtree(temp_dir)
					os.mkdir(temp_dir)
				return outdir, temp_dir, checkpoint_temp
			if os.path.isdir(outdir) and os.listdir(outdir) != []: #outdir exists, not empty, no force overwrite.
				print "Previous results exist! Run with -f if you wish to overwrite previous results, of -c if you wish to run using checkpoints."
				sys.exit(1)
			elif os.path.isdir(outdir) and os.listdir(outdir) == []: #outdir exists, is empty.
				if os.path.isdir(temp_dir) == True:
					shutil.rmtree(temp_dir)
				os.mkdir(temp_dir)
				return outdir, temp_dir, checkpoint_temp
			else:
				os.mkdir(outdir)
		elif sys.argv[2] == 'true':
			print "Beginning bact_vcp..."
			print "Overwriting any previous results."
			if os.path.isdir(outdir):
				shutil.rmtree(outdir)
			os.mkdir(outdir)
		if os.path.isdir(temp_dir) == True:
			shutil.rmtree(temp_dir)
		os.mkdir(temp_dir)
		return outdir, temp_dir, checkpoint_temp
	except:
		shutil.rmtree(outdir)
		os.rename(temp_dir,checkpoint_temp)
		raise


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
	run_piped_shell_process_outerr("java -jar %s build -genbank -v reference" % snpeff_path,outfile)


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
	try:
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
			print "Indexing referencing sequence for picard-tools:\n"
			run_piped_shell_process_outerr("java -Xmx2g -XX:+UseSerialGC -jar %s CreateSequenceDictionary REFERENCE=reference.fa OUTPUT=reference.dict" % picard_path, outfile)
			print "Indexing referencing sequence for smalt:\n"
			run_piped_shell_process_outerr("smalt index -k 13 -s 8 reference reference.fa" % outfile)
			print "Indexing referencing sequence for gatk:\n"
			run_piped_shell_process_outerr("samtools faidx reference.fa" % outfile)
		shutil.copy(os.path.join(temp_dir,"reference.fa"),os.path.join(outdir))
	except:
		shutil.rmtree(outdir)
		shutil.rmtree(temp_dir)
		raise



def trim_reads(read_data,base_path,temp_dir,outdir):
	trimo_path = os.path.join(base_path,"bin/trimmomatic/trimmomatic-0.36.jar")
	cut_adapt_path = os.path.join(base_path,".venv/bin/cutadapt")
	reads, rgid, sample_name, rgpl, rglb = read_data
	if os.path.isfile(reads) == False:
		print "Error: Supplied reads file does not exist!"
		sys.exit(1)
	results_dir = os.path.join(outdir,sample_name)
	os.mkdir(results_dir)
	paired_output1 = "%s.trimmed_paired1.fastq" % sample_name
	paired_output2 = "%s.trimmed_paired2.fastq" % sample_name
	unpaired_output1 = "%s.trimmed_unpaired1.fastq" % sample_name
	unpaired_output2 = "%s.trimmed_unpaired2.fastq" % sample_name
	merged_trimmed = "%s.trimmed.fastq.gz" % sample_name
	print "%s: Importing reads" % sample_name 
	trim_log = "%s.trimlog" % sample_name
	reads1 = "%s.1.fastq" % sample_name
	reads2 = "%s.2.fastq" % sample_name
	ca_trim_log = "%s.ca_trimlog" % sample_name
	with cd(temp_dir):
		print "%s: Spliting reads file" % sample_name 
		unmerge(reads, sample_name, gz=False)
		print "%s: Trimming reads" % sample_name 
		cut_adapt_out1 = "ca_%s" % reads1
		cut_adapt_out2 = "ca_%s" % reads2
		run_piped_shell_process("%s -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o %s -p %s %s %s" % (cut_adapt_path, cut_adapt_out1, cut_adapt_out2, reads1, reads2), ca_trim_log)
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",trimo_path,"PE","-phred33","-trimlog",trim_log,cut_adapt_out1,cut_adapt_out2,paired_output1,unpaired_output1,paired_output2,unpaired_output2,"LEADING:20","TRAILING:20","SLIDINGWINDOW:4:15","MINLEN:70"])
		print "%s: Merging trimmed reads" % sample_name
		with gzip.open(merged_trimmed, 'wb') as trimmed_merged:
			fastqs = [FASTQ(x) for x in (paired_output1,paired_output2)]
			merge(fastqs, split_slashes=False, out=trimmed_merged, quiet=False)
		os.remove(reads1)
		os.remove(reads2)
		os.remove(cut_adapt_out1)
		os.remove(cut_adapt_out2)
	shutil.copyfile(os.path.join(temp_dir,merged_trimmed),os.path.join(results_dir,merged_trimmed))


def sam_to_bam(sam,read_data,base_path):
	reads, rgid, sample_name, rgpl, rglb = read_data
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	bam = "%s.raw_alignment.bam" % sample_name
	pre_bam = "%s.raw_no_rg.bam" % sample_name
	run_piped_shell_process("samtools view -bh -q 10 %s" % sam, pre_bam)
	run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"AddOrReplaceReadGroups","I=%s" % pre_bam,"O=%s" % bam,"RGID=%s" % rgid,"RGSM=%s" % sample_name,"RGPL=%s" % rgpl.lower(),"RGLB=%s" % rglb,"RGPU=NA","SORT_ORDER=coordinate"])
	os.remove(sam)
	os.remove(pre_bam)
	alignment_quality_check(bam,picard_path,sample_name)
	return bam


def alignment_quality_check(bam,picard_path,sample_name):
	stats_dir = "alignment_quality_stats"
	os.mkdir(stats_dir)
	prefix = os.path.join(os.path.abspath(stats_dir),sample_name)
	run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"CollectGcBiasMetrics","R=reference.fa", "I=%s" % bam, "O=%s_GCBias.txt" % prefix, "CHART=%s_GCBias.pdf" % prefix, "ASSUME_SORTED=true"])
	run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"MeanQualityByCycle","R=reference.fa", "I=%s" % bam, "O=%s_Qcycle.txt" % prefix, "CHART=%s_Qcycle.pdf" % prefix])
	run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"QualityScoreDistribution","R=reference.fa", "I=%s" % bam, "O=%s_Qdist.txt" % prefix, "CHART=%s_Qdist.pdf" % prefix])
	




def run_smalt(read_data,base_path,temp_dir,outdir,merged_trimmed):
	reads, rgid, sample_name, rgpl, rglb = read_data
	results_dir = os.path.join(outdir,sample_name)
	reads1 = "%s.1.fastq" % sample_name
	reads2 = "%s.2.fastq" % sample_name
	with cd(temp_dir):
		if os.path.isfile(reads1) == False or os.path.isfile(reads2) == False:
			print "Reads file not found in temp directory. Copying from results directory."
			if os.path.isfile(reads1) == True:
				os.remove(reads1)
			if os.path.isfile(reads2) == True:
				os.remove(reads2)
			unmerge(merged_trimmed, sample_name, gz=False)
		print "%s: Processing reads for mapping" % sample_name
		sam = "%s.sam" % sample_name
		unmerge(reads, sample_name, gz=False)
		print "%s: Mapping reads" % sample_name
		run_process(["smalt","map","-i","1500","-n","1","-o",sam,"reference",reads1,reads2])
		print "%s: Generating raw alignment bam file" % sample_name
		bam = sam_to_bam(sam,read_data,base_path)
	shutil.move(os.path.join(temp_dir,"alignment_quality_stats"),os.path.join(results_dir))
	shutil.copyfile(os.path.join(temp_dir,bam),os.path.join(results_dir,bam))



#Counts number of snps in a 200 bp window across whole genome
def count_variants_and_plot(filtered_variants):
	count_list = []
	win_out = []
	window = 1000
	records = [record for record in vcf.Reader(open(filtered_variants, 'r'))]
	last_entry = records[-1].POS
	snp_positions = numpy.zeros(last_entry+1, dtype=bool)
	for record in records:
		snp_positions[record.POS] = True
	temp_counts = []
	for n in xrange(0, last_entry, window):
		count = numpy.count_nonzero(snp_positions[n:n+window])
		temp_counts.append(count)
		if i == 0:
			win_out.append(n)
	if i == 0:
		count_list.append(win_out)	
	count_list.append(temp_counts)

	with open("snp_density.csv", "wb") as fh_out:
		writer = csv.writer(fh_out)
		files.insert(0, 'position')
		writer.writerow(files)
		writer.writerows(zip(*count_list))
	


def dedup_reads(read_data, base_path,temp_dir,outdir,results_bam):
	print "%s: Removing PCR duplicates" % read_data[2]
	reads, rgid, sample_name, rgpl, rglb = read_data
	results_dir = os.path.join(outdir,sample_name)
	bam = os.path.join(temp_dir,"%s.raw_alignment.bam" % sample_name)
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	dedup_bam = "%s.dedup_reads.bam" % sample_name
	markd_input = "INPUT=%s" % bam
	markd_output = "OUTPUT=%s" % dedup_bam
	markd_metrics = "%s.picard_md_metrics.txt" % sample_name
	markd_metrics_arg = "METRICS_FILE=%s" % markd_metrics
	bbiinput = "INPUT=%s" % dedup_bam
	with cd(temp_dir):
		if os.path.isfile(bam) == False:
			print "Bam alignment file not found in temp directory. Copying from results directory."
			shutil.copyfile(results_bam,bam)
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"MarkDuplicates",markd_input,markd_output,markd_metrics_arg])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput])
	shutil.copyfile(os.path.join(temp_dir,dedup_bam),os.path.join(results_dir,dedup_bam))
	shutil.copyfile(os.path.join(temp_dir,markd_metrics),os.path.join(results_dir,markd_metrics))


def realign_indels(read_data,base_path,temp_dir,outdir,results_dedup_bam):
	print "%s: Realigning indel mutations" % read_data[2]
	reads, rgid, sample_name, rgpl, rglb = read_data
	results_dir = os.path.join(outdir,sample_name)
	dedup_bam = "%s.dedup_reads.bam" % sample_name
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	realign_intervals = "%s.target.intervals" % sample_name
	preprocessed_reads = "%s.preprocessed_reads.bam" % sample_name
	bbiinput = "INPUT=%s" % preprocessed_reads
	cov_stats_name = "%s.coverage_stats" % sample_name
	with cd(temp_dir):
		if os.path.isfile(dedup_bam) == False:
			print "Deduplicated bam alignment file not found in temp directory. Copying from results directory."
			shutil.copyfile(results_dedup_bam,dedup_bam)
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","RealignerTargetCreator","-R","reference.fa","-I",dedup_bam,"-o",realign_intervals])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","IndelRealigner","-R","reference.fa","-I",dedup_bam,"-targetIntervals",realign_intervals,"-o",preprocessed_reads])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput])
		print "%s: Calculating depth of coverage statistics" % sample_name
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","DepthOfCoverage","-R","reference.fa","-o",cov_stats_name,"-I",preprocessed_reads,"-ct","5","-ct","10"])
	shutil.copyfile(os.path.join(temp_dir,preprocessed_reads),os.path.join(results_dir,preprocessed_reads))


def call_variants(read_data,base_path,temp_dir,outdir,results_preprocessed_reads):
	print "%s: Calling variants" % read_data[2]
	reads, rgid, sample_name, rgpl, rglb = read_data
	results_dir = os.path.join(outdir,sample_name)
	preprocessed_reads = "%s.preprocessed_reads.bam" % sample_name
	gatk_path = os.path.join(base_path,"gatk/GenomeAnalysisTK.jar")
	raw_variants = "%s.raw_variants.vcf" % sample_name
	with cd(temp_dir):
		if os.path.isfile(preprocessed_reads) == False:
			print "Preprocessed bam alignment file not found in temp directory. Copying from results directory."
			shutil.copyfile(results_preprocessed_reads,preprocessed_reads)
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","HaplotypeCaller","-R","reference.fa","-I",preprocessed_reads,"--genotyping_mode","DISCOVERY","-stand_call_conf","30","-stand_emit_conf","10","-ERC", "GVCF", "-o",raw_variants,"--sample_ploidy","1"])
	shutil.copyfile(os.path.join(temp_dir,raw_variants),os.path.join(results_dir,raw_variants))


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
	with cd(temp_dir):
		if os.path.isfile(raw_variants) == False:
			print "Raw variants file not found in temp directory. Copying from results directory."
			shutil.copyfile(os.path.join(outdir,raw_variants),raw_variants)
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",raw_variants,"-selectType","SNP","-o",raw_snps])
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",raw_variants,"-selectType","INDEL","-o",raw_indels])
		cmd = 'java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s ' % (gatk_path,raw_snps)
		cmd += '--filterExpression "QD < 2.0" --filterName "low_qual_by_depth" '
		cmd += '--filterExpression "MQ < 40.0" --filterName "low_RMS_mapping_quality" '
		cmd += '--filterExpression "FS > 60.0" --filterName "strand_bias" '
		cmd += '--filterExpression "SOR > 3.0" --filterName "high_strand_odds_ratio" '
		cmd += '--filterExpression "MQRankSum < -12.5" --filterName "low_MQRankSum" '
		cmd += '--filterExpression "ReadPosRankSum < -8.0" --filterName "low_ReadPosRankSum" -o %s' % filtered_snps
		run_shell_process(cmd)
		cmd = 'java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s ' % (gatk_path,raw_indels)
		cmd += '--filterExpression "QD < 2.0" --filterName "low_qual_by_depth" '
		cmd += '--filterExpression "ReadPosRankSum < -20.0" --filterName "low_ReadPosRankSum" '
		cmd += '--filterExpression "FS > 200.0" --filterName "strand_bias" '
		cmd += '--filterExpression "SOR > 10.0" --filterName "high_strand_odds_ratio" -o %s' % filtered_indels
		run_shell_process(cmd)
		print "%s: Annotating variants" % read_data[2]
		run_process(["java", "-Xmx2g", "-XX:+UseSerialGC", "-jar", gatk_path, "-T", "CombineVariants", "-R", "reference.fa", "--variant:snp", filtered_snps, "--variant:indel", filtered_indels, "-o", filtered_variants, "-genotypeMergeOptions", "PRIORITIZE", "-priority", "snp,indel"])
		run_piped_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -v reference %s' % (snpeff_path,filtered_variants),annotated_variants)
	shutil.copyfile(os.path.join(temp_dir,filtered_variants),os.path.join(outdir,filtered_variants))
	shutil.copyfile(os.path.join(temp_dir,annotated_variants),os.path.join(outdir,annotated_variants))





def run_pipeline(read_data,base_path,temp_dir,outdir):
	filtered_variants_file = []
	results_dir = os.path.join(outdir,read_data[2])
	merged_trimmed = os.path.join(results_dir,"%s.trimmed.fastq.gz" % read_data[2])
	bam = os.path.join(results_dir,"%s.raw_alignment.bam" % read_data[2])
	dedup_bam = os.path.join(results_dir,"%s.dedup_reads.bam" % read_data[2])
	preprocessed_reads = os.path.join(results_dir,"%s.preprocessed_reads.bam" % read_data[2])
	raw_variants = os.path.join(results_dir,"%s.raw_variants.vcf" % read_data[2])
	if sys.argv[6] == 'true':
		if os.path.exists(merged_trimmed) == False:
			trim_reads(read_data,base_path,temp_dir,outdir)
		if os.path.exists(bam) == False:
			run_smalt(read_data,base_path,temp_dir,outdir,os.path.abspath(merged_trimmed))
		if os.path.exists(dedup_bam) == False:
			dedup_reads(read_data, base_path,temp_dir,outdir,os.path.abspath(bam))
		if os.path.exists(preprocessed_reads) == False:
			realign_indels(read_data,base_path,temp_dir,outdir,os.path.abspath(dedup_bam))
		if os.path.exists(raw_variants) == False:
			call_variants(read_data,base_path,temp_dir,outdir,os.path.abspath(preprocessed_reads))
	else:
		trim_reads(read_data,base_path,temp_dir,outdir)
		run_smalt(read_data,base_path,temp_dir,outdir,os.path.abspath(merged_trimmed))
		dedup_reads(read_data, base_path,temp_dir,outdir,os.path.abspath(bam))
		realign_indels(read_data,base_path,temp_dir,outdir,os.path.abspath(dedup_bam))
		call_variants(read_data,base_path,temp_dir,outdir,os.path.abspath(preprocessed_reads))


def multiprocessing(read_data):
	outdir = sys.argv[5]
	base_path = os.path.dirname(os.path.dirname(__file__))
	temp_dir = os.path.join(base_path,"intermediate_files")
	run_pipeline(read_data,base_path,temp_dir,outdir)

def genotype_variants(reads_list,temp_dir,outdir):
	variant_files = [os.path.join(outdir,x[2],"%s.raw_variants.vcf" % x[2] for x in reads_list]
	cmd = ["java", "-Xmx2g", "-XX:+UseSerialGC", "-jar", gatk_path, "-T", "GenotypeGVCFs", "-R", "reference.fa"]
	for raw_variants in variant_files:
		cmd.append("--variant")
		cmd.append(raw_variants)
	cmd.append("-o")
	cmd.append("raw_variants.vcf")
	with cd(temp_dir):
		run_process(cmd)
	shutil.copy(os.path.join(temp_dir,"all.raw_variants.vcf"),outdir)


base_path = os.path.dirname(os.path.dirname(__file__))
file_list = sys.argv[1]
ref_gb = sys.argv[3]
outdir = sys.argv[5]



if __name__ == "__main__":
	outdir, temp_dir, checkpoint_temp = setup_directories_and_inputs(base_path, outdir)
	reads_list = import_file_list(file_list,outdir)
	null_output = process_reference(base_path,ref_gb,temp_dir,outdir)
	setup_snpEff(ref_gb,base_path,outdir)
	try:
		#for read_data in reads_list:
			#multiprocessing(read_data)
		p = Pool(int(sys.argv[7]))
		p.map_async(multiprocessing,reads_list).get(9999999)
		p.close()
		p.join()
		genotype_variants(reads_list,temp_dir,outdir)
		annotate_variants(base_path,temp_dir,outdir)
	except:
		if sys.argv[6] == 'true':
			print "Error!!! Saving temp files to checkpoint directory..."
			shutil.move(temp_dir,checkpoint_temp)
		else:
			if sys.argv[4] == 'true':
				print "Error!!! Saving temp files to output directory as requested..."
				shutil.move(temp_dir,os.path.join(outdir,"temp_files"))
			else:
				print "Error!!! Deleting temp files as requested..."
				shutil.rmtree(temp_dir)
		reset_snpEff(ref_gb,base_path)
		raise
	else:
		if sys.argv[4] == 'true':
			print "Run completed successfully. Saving temp files to output directory..."
			if os.path.isdir(os.path.join(outdir,"temp_files")) == True:
				print "Removing temp_files from previous successful run"
				shutil.rmtree(os.path.join(outdir,"temp_files"))
			shutil.move(temp_dir,os.path.join(outdir,"temp_files"))
		elif sys.argv[4] == 'false':
			print "Run completed successfully. Deleting temp files..."
			shutil.rmtree(temp_dir)
		reset_snpEff(ref_gb,base_path)
		print "All jobs completed successfully."


