
VENV = .venv
export VIRTUAL_ENV := $(abspath ${VENV})
export PKG_CONFIG_PATH := ${VIRTUAL_ENV}/lib/pkgconfig
export PATH := ${VIRTUAL_ENV}/bin:${PATH}


all: bin/gatk/GenomeAnalysisTK.jar ${VENV}/bin/python \
${VENV}/bin/smalt \
${VENV}/bin/samtools \
${VENV}/lib/libz.so \
${VENV}/lib/libpng.so \
${VENV}/lib/libfreetype.so \
${VENV}/bin/pkg-config \
${VENV}/bin/cutadapt \

#Add make rule for gatk, use a variable for GATK location

clean: 
	rm -rf bin/smalt-0.7.6 bin/samtools-1.3.1 bin/gatk ${VENV} intermediate_files checkpoint_files

.PHONY: all clean
.SECONDARY:

bin/gatk/GenomeAnalysisTK.jar:
	mkdir bin/gatk
	cp -av ${GATK_PATH}/GenomeAnalysisTK.jar bin/gatk/

${VENV}/bin/python:
	python2 support/virtualenv.py ${VENV}

${VENV}/bin/smalt: ${VENV}/bin/python
	cd bin && tar -xvf smalt-0.7.6-static.tar.gz
	cd bin/smalt-0.7.6 && ./configure --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install

${VENV}/bin/samtools: ${VENV}/bin/python
	cd bin && tar -xvf samtools-1.3.1.tar.bz2
	cd bin/samtools-1.3.1 && ./configure --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install

${VENV}/lib/libz.so: ${VENV}/bin/python
	cd bin && tar -xvf zlib-1.2.8.tar.gz
	cd bin/zlib-1.2.8 && ./configure --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install
	rm -rf bin/zlib-1.2.8

${VENV}/lib/libpng.so: ${VENV}/bin/python
	cd bin && tar -xvf libpng-1.6.24.tar.gz
	cd bin/libpng-1.6.24 && ./configure --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install
	rm -rf bin/libpng-1.6.24

${VENV}/lib/libfreetype.so: ${VENV}/bin/python
	cd bin && tar -xvf freetype-2.6.5.tar.gz
	cd bin/freetype-2.6.5 && ./configure --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install
	rm -rf bin/freetype-2.6.5

${VENV}/bin/pkg-config: ${VENV}/bin/python
	cd bin && tar -xvf pkg-config-0.29.tar.gz
	cd bin/pkg-config-0.29 && ./configure --prefix=${VIRTUAL_ENV} --with-internal-glib && $(MAKE) && $(MAKE) install
	rm -rf bin/pkg-config-0.29


${VENV}/bin/cutadapt: ${VENV}/bin/python
	pip install --upgrade -r scripts/requirements.pip
