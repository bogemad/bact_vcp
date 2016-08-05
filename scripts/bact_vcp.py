#!/usr/bin/env python

import os, sys, shutil, subprocess, collections, gzip, re
from eta import ETA
from contextlib import contextmanager


base_path = os.path.dirname(os.path.dirname(__file__))

file_list = sys.argv[1]
force = sys.argv[2]
ref_fasta = sys.argv[3]
ref_db = sys.argv[4]
keep = sys.argv[5]
outdir = sys.argv[6]

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
        if self.fname and not quiet:
            eta = ETA(os.stat(self.fname).st_size, fileobj=self.fileobj)
        else:
            eta = None

        while True:
            try:
                read = fastq_read_file(self.fileobj)
                if eta:
                    if callback:
                        eta.print_status(extra=callback())
                    else:
                        eta.print_status(extra=read.name)
                yield read

            except StopIteration:
                break

        if eta:
            eta.done()

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
	output = subprocess.check_output(cmd, stderr = subprocess.STDOUT)
	print output


def run_shell_process(cmd):
	output = subprocess.check_output(cmd, stderr = subprocess.STDOUT, shell=True)
	print output


def run_piped_shell_process(cmd):
	subprocess.call(cmd, shell=True)


def setup_directories_and_inputs(base_path, outdir, force):
	temp_dir = os.path.join(base_path,"intermediate_files")
	if force == 'false':
		if os.path.isdir(outdir) and os.listdir(outdir) != []: #outdir exists, not empty, no force overwrite.
			print "Previous results exist! Run with -f if you wish to overwrite previous results."
			sys.exit(1)
		elif os.path.isdir(outdir) and os.listdir(outdir) == []: #outdir exists, is empty.
			
			if os.path.isdir(temp_dir) == True:
				shutil.rmtree(temp_dir)
			os.mkdir(temp_dir)
			return outdir, temp_dir
		else:
			os.mkdir(outdir)
	elif force == 'true':
		print "Beginning bact_vcp..."
		print "Overwriting any previous results."
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		os.mkdir(outdir)
	if os.path.isdir(temp_dir) == True:
		shutil.rmtree(temp_dir)
	os.mkdir(temp_dir)
	return outdir, temp_dir


def import_file_list(file_list):
	reads_list = []
	with open(file_list) as file_handle:
		for line in file_handle:
			line = re.sub(r'^"#.*"',"",line) # fixes quotes that MS excel sometimes adds to csv files
			if line.startswith("#") or line.startswith(",") or line == "\n":
				continue
			data = line.rstrip().split(',')
			reads_list.append((data[0],data[1],data[2],data[3],data[4]))
	return reads_list


def copy_reference(base_path,ref_fasta,temp_dir,outdir):
	picard_path = os.path.join(base_path, "bin/picard-tools/picard.jar")
	print "Importing reference sequence: %s" % os.path.basename(ref_fasta)
	ref_path = os.path.join(temp_dir,"reference.fa")
	shutil.copyfile(ref_fasta,ref_path)
	with cd(temp_dir):
		if os.path.isfile("reference.dict"): os.remove("reference.dict")
		run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"CreateSequenceDictionary","REFERENCE=reference.fa","OUTPUT=reference.dict"])
	shutil.copy(os.path.join(temp_dir,"reference.fa"),os.path.join(outdir))


def trim_reads(reads_list, base_path,temp_dir,outdir):
	trimo_path = os.path.join(base_path,"bin/trimmomatic/trimmomatic-0.36.jar")
	cut_adapt_path = os.path.join(base_path,".venv/bin/cutadapt")
	for read_data in reads_list:
		reads, rgid, sample_name, rgpl, rglb = read_data
		results_dir = os.path.join(outdir,sample_name)
		paired_output1 = "%s.trimmed_paired1.fastq" % sample_name
		paired_output2 = "%s.trimmed_paired2.fastq" % sample_name
		unpaired_output1 = "%s.trimmed_unpaired1.fastq" % sample_name
		unpaired_output2 = "%s.trimmed_unpaired2.fastq" % sample_name
		merged_trimmed = "%s.trimmed.fastq.gz" % sample_name
		print "%s: Importing reads" % sample_name 
		trim_log = "%s.trimlog" % sample_name
		ca_trim_log = "%s.ca_trimlog" % sample_name
		reads1 = "%s.1.fastq" % sample_name
		reads2 = "%s.2.fastq" % sample_name
		with cd(temp_dir):
			print "%s: Spliting reads file" % sample_name 
			unmerge(reads, sample_name, gz=False)
			print "%s: Trimming reads" % sample_name 
			cut_adapt_out1 = "ca_%s" % reads1
			cut_adapt_out2 = "ca_%s" % reads2
			run_piped_shell_process("%s -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o %s -p %s %s %s > %s" % (cut_adapt_path, cut_adapt_out1, cut_adapt_out2, reads1, reads2, ca_trim_log))
			run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",trimo_path,"PE","-phred33","-trimlog",trim_log,cut_adapt_out1,cut_adapt_out2,paired_output1,unpaired_output1,paired_output2,unpaired_output2,"LEADING:20","TRAILING:20","SLIDINGWINDOW:4:15","MINLEN:70"])
			with gzip.open(os.path.join(results_dir,merged_trimmed), 'w') as trimmed_merged:
				fastqs = [FASTQ(x) for x in (paired_output1,paired_ooutput2)]
				merge(fastqs, split_slashes=False, out=trimmed_merged, quiet=False)
			os.remove(reads1)
			os.remove(reads2)
			os.remove(cut_adapt_out1)
			os.remove(cut_adapt_out2)


def sam_to_bam(sam,sample_name,base_path):
	picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
	bam = "%s.raw_alignment.bam" % sample_name
	bam_arg = "O=%s" % bam
	samfile_arg = "I=%s" % sam
	run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"SortSam",samfile_arg,bam_arg,"SORT_ORDER=coordinate"])
	os.remove(sam)


#read-groups!
def run_smalt(reads_list,base_path,temp_dir):
	for read_data in reads_list:
		reads, rgid, sample_name, rgpl, rglb = read_data
		results_dir = os.path.join(outdir,sample_name)
		trimmed_reads = os.path.join(results_dir,"%s.trimmed.fastq.gz" % sample_name)
		reads1 = "%s.1.fastq" % sample_name
		reads2 = "%s.2.fastq" % sample_name
		smalt_path = os.path.join(base_path,"bin/smalt-0.7.6/src/smalt")
		with cd(temp_dir):
			print "%s: Mapping reads" % sample_name
			run_process([smalt_path,"index","-k","13","-s","8","reference","reference.fa"])
			sam = "%s.sam" % sample_name
			unmerge(reads, sample_name, gz=False)
			run_process([smalt_path,"map","-i","1500","-n","1","-o",sam,"reference",reads1,reads2])
			print "%s: Generating bam output" % sample_name
			bam = sam_to_bam(sam,sample_name,base_path)
			shutil.copyfile(bam,os.path.join(results_dir,bam))


def dedup_reads(reads_list, base_path, temp_dir,outdir):
	for read_data in reads_list:
		print "%s: Removing PCR duplicates" % read_data[2]
		reads, rgid, sample_name, rgpl, rglb = read_data
		results_dir = os.path.join(outdir,sample_name)
		bam = os.path.join(temp_dir,"%s.raw_alignment.bam" % sample_name)
		picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
		dedup_bam = "%s.dedup_reads.bam" % sample_name
		markd_input = "INPUT=%s" % bam
		markd_output = "OUTPUT=%s" % dedup_bam
		markd_metrics = "METRICS_FILE=%s.picard_md_metrics.txt" % os.path.join(outdir,sample_name,sample_name)
		bbiinput = "INPUT=%s" % dedup_bam
		with cd(temp_dir):
			run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"MarkDuplicates",markd_input,markd_output,markd_metrics])
			run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput])
			shutil.copyfile(dedup_bam, os.path.join(results_dir,dedup_bam))


def realign_indels(reads_list, base_path,temp_dir,outdir):
	for read_data in reads_list:
		print "%s: Realigning indel mutations" % read_data[2]
		reads, rgid, sample_name, rgpl, rglb = read_data
		results_dir = os.path.join(outdir,sample_name)
		dedup_bam = "%s.dedup_reads.bam" % sample_name
		picard_path = os.path.join(base_path,"bin/picard-tools/picard.jar")
		gatk_path = os.path.join(base_path,"bin/gatk/GenomeAnalysisTK.jar")
		realign_intervals = "%s.target_intervals.txt" % sample_name
		preprocessed_reads = "%s.preprocessed_reads.bam" % sample_name
		bbiinput = "INPUT=%s" % preprocessed_reads
		with cd(temp_dir):
			run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","RealignerTargetCreator","-R","reference.fa","-I",dedup_bam,"-o",realign_intervals])
			run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","IndelRealigner","-R","reference.fa","-I",dedup_bam,"-targetIntervals",realign_intervals,"-o",preprocessed_reads])
			run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",picard_path,"BuildBamIndex",bbiinput])
			shutil.copyfile(preprocessed_reads, os.path.join(results_dir,preprocessed_reads))


def call_variants(reads_list, base_path,temp_dir,outdir):
	for read_data in reads_list:
		print "%s: Calling variants" % read_data[2]
		reads, rgid, sample_name, rgpl, rglb = read_data
		results_dir = os.path.join(outdir,sample_name)
		preprocessed_reads = "%s.preprocessed_reads.bam" % sample_name
		gatk_path = os.path.join(base_path,"bin/gatk/GenomeAnalysisTK.jar")
		raw_variants = "%s.raw_variants.vcf" % sample_name
		with cd(temp_dir):
			run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","HaplotypeCaller","-R","reference.fa","-I",preprocessed_reads,"--genotyping_mode","DISCOVERY","-glm","BOTH","-stand_call_conf","30","-stand_emit_conf","10","-o",raw_variants,"--sample_ploidy","1"])
			shutil.copyfile(raw_variants, os.path.join(results_dir,raw_variants))


def annotate_snps(reads_list,base_path,temp_dir,outdir,ref_db):
	for read_data in reads_list:
		print "%s: First step filtering of SNPs" % read_data[2]
		reads, rgid, sample_name, rgpl, rglb = read_data
		results_dir = os.path.join(outdir,sample_name)
		raw_variants = "%s.raw_variants.vcf" % sample_name
		gatk_path = os.path.join(base_path,"bin/gatk/GenomeAnalysisTK.jar")
		raw_snps = "%s.raw_snps.vcf" % sample_name
		filtered_snps = "%s.filtered_snps.vcf" % sample_name
		snpeff_path = os.path.join(base_path,"bin/snpEff/snpEff.jar")
		annotated_snps = "%s.annotated_snps.vcf" % sample_name
		with cd(temp_dir):
			run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",raw_variants,"-selectType","SNP","-o",raw_snps])
			run_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s --filterExpression "QD < 3.0" --filterName "low_qual_by_depth" --filterExpression "FS > 13.5" --filterName "strand_bias" --filterExpression "DP < 10" --filterName "low_depth" --filterExpression "MQRankSum < -12.5" --filterName "MQRankSum" --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum" -o %s' % (gatk_path,raw_snps,filtered_snps))
			run_piped_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -v %s %s > %s' % (snpeff_path,ref_db,filtered_snps,annotated_snps))
			shutil.copyfile(filtered_snps,os.path.join(results_dir,filtered_snps))
			shutil.copyfile(annotated_snps,os.path.join(results_dir,annotated_snps))


def annotate_indels(reads_list,base_path,temp_dir,outdir,ref_db):
	for read_data in reads_list:
		print "%s: First step filtering of Indels" % read_data[2]
		reads, rgid, sample_name, rgpl, rglb = read_data
		results_dir = os.path.join(outdir,sample_name)
		raw_variants = "%s.raw_variants.vcf" % sample_name
		gatk_path = os.path.join(base_path,"bin/gatk/GenomeAnalysisTK.jar")
		raw_indels = "%s.raw_indels.vcf" % sample_name
		filtered_indels = "%s.filtered_indels.vcf" % sample_name
		snpeff_path = os.path.join(script_path,"bin/snpEff/snpEff.jar")
		annotated_indels = "%s.annotated_indels.vcf" % sample_name
		with cd(temp_dir):
			run_process(["java","-Xmx2g","-XX:+UseSerialGC","-jar",gatk_path,"-T","SelectVariants","-R","reference.fa","-V",raw_variants,"-selectType","SNP","-o",raw_indels])
			run_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -T VariantFiltration -R reference.fa -V %s --filterExpression "QD < 3.0" --filterName "low_qual_by_depth" --filterExpression "FS > 200.0" --filterName "strand_bias" --filterExpression "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum" --filterExpression "DP < 10" --filterName "low_depth" -o %s' % (gatk_path,raw_indels,filtered_indels))
			run_piped_shell_process('java -Xmx2g -XX:+UseSerialGC -jar %s -v %s %s > %s' % (snpeff_path,ref_db,filtered_indels,annotated_indels))
			shutil.copyfile(filtered_indels,os.path.join(results_dir,filtered_indels))
			shutil.copyfile(annotated_snps,os.path.join(results_dir,annotated_snps))






if __name__ == "__main__":
	outdir, temp_dir = setup_directories_and_inputs(base_path, outdir, force)
	reads_list = import_file_list(file_list)
	null_output = copy_reference(base_path,ref_fasta,temp_dir,outdir)
	trim_reads(reads_list, base_path,temp_dir,outdir)
	run_smalt(reads_list,base_path,temp_dir,sample_name)
	dedup_reads(reads_list, base_path, temp_dir,outdir)
	realign_indels(reads_list, base_path,temp_dir,outdir)
	call_variants(reads_list, base_path,temp_dir,outdir)
	annotate_snps(reads_list,base_path,temp_dir,outdir,ref_db)
	annotate_indels(reads_list,base_path,temp_dir,outdir,ref_db)
	if keep == 'true':
		shutil.move(temp_dir,os.path.join(outdir,"temp_files"))
	elif keep == 'false':
		shutil.rmtree(temp_dir)
	else:
		print "%s <-- This is the keep value" % keep
