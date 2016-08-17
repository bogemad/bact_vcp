
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


all: gatk/GenomeAnalysisTK.jar ${VENV}/bin/python \
${VENV}/bin/smalt \
${VENV}/lib/libz.so \
${VENV}/bin/samtools \

#Add make rule for gatk, use a variable for GATK location

clean: 
	rm -rf smalt-0.7.6-static.tar.gz smalt-0.7.6 zlib-1.2.8 zlib-1.2.8.tar.gz samtools-1.3.1 samtools-1.3.1.tar.bz2 gatk mc.sh ${MC} intermediate_files checkpoint_files

.PHONY: all clean
.SECONDARY:

gatk/GenomeAnalysisTK.jar:
	mkdir gatk
	cp -av ${GATK_PATH}/GenomeAnalysisTK.jar gatk/

${VENV}/bin/python:
	wget -O mc.sh ${MC_LINK}
	chmod 755 mc.sh && ./mc.sh -b -p ${MC}
	${MC}/bin/conda create -y -n venv python=2.7.11 pip numpy matplotlib biopython
	${VENV}/bin/pip install pyvcf cutadapt
	rm -f mc.sh

${VENV}/bin/smalt: ${VENV}/bin/python
	wget http://downloads.sourceforge.net/project/smalt/smalt-0.7.6-static.tar.gz 
	tar -xvf smalt-0.7.6-static.tar.gz
	cd smalt-0.7.6 && ./configure --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install
	rm -rf smalt-0.7.6-static.tar.gz smalt-0.7.6 

${VENV}/lib/libz.so: ${VENV}/bin/python
	wget http://downloads.sourceforge.net/project/libpng/zlib/1.2.8/zlib-1.2.8.tar.gz
	tar -xvf zlib-1.2.8.tar.gz
	cd zlib-1.2.8 && ./configure --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install
	rm -rf zlib-1.2.8 zlib-1.2.8.tar.gz

${VENV}/bin/samtools: ${VENV}/bin/python ${VENV}/lib/libz.so
	wget -O samtools-1.3.1.tar.bz2 https://sourceforge.net/projects/samtools/files/samtools/1.3.1/samtools-1.3.1.tar.bz2/download
	tar -xvf samtools-1.3.1.tar.bz2
	cd samtools-1.3.1 && ./configure --without-curses --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install
	rm -rf samtools-1.3.1 samtools-1.3.1.tar.bz2
