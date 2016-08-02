all: bin/smalt-0.7.6/src/smalt bin/ngsutils-ngsutils-0.5.9/venv/bin/activate ${VENV} python-reqs make_keyfiles ${RESULTS}/%.trimmed.fastq.gz ${RESULTS}/%.raw_alignment.bam ${RESULTS}/%.dedup_reads.bam ${RESULTS}/%.gatk_preprocessed_reads.bam ${RESULTS}/%.raw_variants.vcf ${RESULTS}/%.annotated_snps.vcf ${RESULTS}/%.annotated_indels.vcf

export KEYS := $(abspath keys)
export RESULTS := $(abspath results)
export SCRIPTS := $(abspath scripts)
export INT_FILES := $(abspath intermediate_files)

VENV = .venv
export VIRTUAL_ENV := $(abspath ${VENV})
export PATH := ${VIRTUAL_ENV}/bin:${PATH}

clean: 
	rm -rf bin/smalt-0.7.6 bin/ngsutils-ngsutils-0.5.9 ${VENV} ${KEYS} ${RESULTS} ${INT_FILES}

.PHONY: all clean
.SECONDARY:

bin/smalt-0.7.6/src/smalt:
	cd bin && tar -xvf smalt-0.7.6-static.tar.gz
	cd bin/smalt-0.7.6 && ./configure && $(MAKE)

bin/ngsutils-ngsutils-0.5.9/venv/bin/activate:
	cd bin && tar -xvf ngsutils-0.5.9.tar.gz
	cd bin/ngsutils-ngsutils-0.5.9 && $(MAKE)

${VENV}: bin/smalt-0.7.6/src/smalt bin/ngsutils-ngsutils-0.5.9/venv/bin/activate
	python2 -m virtualenv $@

python-reqs: scripts/requirements.pip | ${VENV}
	pip install --upgrade -r scripts/requirements.pip

${KEYS} ${RESULTS} ${INT_FILES}:
	mkdir $@


make_keyfiles: ${VENV} python-reqs | ${KEYS} ${RESULTS} ${INT_FILES}
	${SCRIPTS}/process_references_lists.py input_data_list.csv 


#maybe make target a compressed archive to simply for below rule
${RESULTS}/%.trimmed.fastq.gz: ${SCRIPTS}/trim_reads.py ${KEYS}/%.txt 
	$^ $@

${RESULTS}/%.raw_alignment.bam: ${SCRIPTS}/gather_data_align.py ${KEYS}/%.txt ${RESULTS}/%.trimmed.fastq.gz
	$^ $@

${RESULTS}/%.dedup_reads.bam: ${SCRIPTS}/dedup_reads.py ${KEYS}/%.txt ${RESULTS}/%.raw_alignment.bam
	$^ $@

${RESULTS}/%.gatk_preprocessed_reads.bam: ${SCRIPTS}/realign_indels.py ${KEYS}/%.txt ${RESULTS}/%.dedup_reads.bam
	$^ $@

${RESULTS}/%.raw_variants.vcf: ${SCRIPTS}/call_variants.py ${KEYS}/%.txt ${RESULTS}/%.gatk_preprocessed_reads.bam
	$^ $@

${RESULTS}/%.annotated_snps.vcf: ${SCRIPTS}/annotate_snps.py ${KEYS}/%.txt ${RESULTS}/%.raw_variants.vcf
	$^ $@

${RESULTS}/%.annotated_indels.vcf: ${SCRIPTS}/annotate_indels.py ${KEYS}/%.txt ${RESULTS}/%.raw_variants.vcf
	$^ $@
