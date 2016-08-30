
MC = .mc
VENV = .mc/envs/venv
export VIRTUAL_ENV := $(abspath ${VENV})
export PATH := ${VIRTUAL_ENV}/bin:${PATH}

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CCFLAGS += -D LINUX
	MC_LINK := https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
endif
ifeq ($(UNAME_S),Darwin)
	CCFLAGS += -D OSX
	MC_LINK := https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
endif
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),x86_64)
	CCFLAGS += -D AMD64
endif
ifneq ($(filter %86,$(UNAME_P)),)
	CCFLAGS += -D IA32
endif
ifneq ($(filter arm%,$(UNAME_P)),)
	CCFLAGS += -D ARM
endif

check_defined = \
    $(strip $(foreach 1,$1, \
        $(call __check_defined,$1,$(strip $(value 2)))))
__check_defined = \
    $(if $(value $1),, \
      $(error Undefined $1$(if $2, ($2))))


all: gatk/GenomeAnalysisTK.jar ${VENV}/bin/python

clean: 
	rm -rf scripts/*.pyc gatk ${MC} intermediate_files checkpoint_files scripts/mc.sh

.PHONY: all clean
.SECONDARY:

gatk/GenomeAnalysisTK.jar:
	$(call check_defined, GATK_PATH, error variable GATK_PATH is not set - set with command: export GATK_PATH=/path/to/gatk_directory)
	mkdir gatk
	cp -av ${GATK_PATH}/GenomeAnalysisTK.jar gatk/
	sed "s~export GATK_PATH=~export GATK_PATH=${GATK_PATH}~g" scripts/hpc.job_lib.sh > scripts/hpc.job_lib2.sh
	rm scripts/hpc.job_lib.sh && mv scripts/hpc.job_lib2.sh scripts/hpc.job_lib.sh

${VENV}/bin/python:
	wget -O - ${MC_LINK} > scripts/mc.sh 
	chmod 755 -R hpc.data.list.csv scripts
	scripts/mc.sh -b -p ${MC}
	${MC}/bin/conda config --add channels r
	${MC}/bin/conda config --add channels bioconda 
	${MC}/bin/conda create -y -n venv python=2.7.11 pip zlib numpy matplotlib biopython psutil smalt samtools java-jdk cutadapt pyvcf
	rm -fr scripts/mc.sh
