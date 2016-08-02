all: ${VENV} python-reqs 

clean: 


.PHONY: all clean
.SECONDARY:

bin/smalt/src/smalt:
	cd bin/smalt && ./configure && $(MAKE)

bin/ngsutils/venv/bin/activate:
	cd bin/ngsutils && $(MAKE)

VENV = .venv
export VIRTUAL_ENV := $(abspath ${VENV})
export PATH := ${VIRTUAL_ENV}/bin:${PATH}


${VENV}: 
	python2 -m virtualenv $@

python-reqs: scripts/requirements.pip | ${VENV}
	pip install --upgrade -r scripts/requirements.pip

export KEYS := $(abspath keys)
export RESULTS := $(abspath results)
export SCRIPTS := $(abspath scripts)

${KEYS}/%.txt: ${SCRIPTS}/process_references_lists.py input_data_list.csv
	mkdir ${KEYS}
	mkdir ${RESULTS}
	mkdir intermediate_files
	$^

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
