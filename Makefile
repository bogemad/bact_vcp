
VENV = .venv
export VIRTUAL_ENV := $(abspath ${VENV})
export PATH := ${VIRTUAL_ENV}/bin:${PATH}
export LDFLAGS := -L${VIRTUAL_ENV}/lib
export CFLAGS := -I${VIRTUAL_ENV}/include
export LD_LIBRARY_PATH := ${VIRTUAL_ENV}/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH := ${VIRTUAL_ENV}/include
export CPLUS_INCLUDE_PATH := ${VIRTUAL_ENV}/include

all: bin/gatk/GenomeAnalysisTK.jar ${VENV}/bin/python \
${VENV}/bin/smalt \
${VENV}/lib/libz.so \
${VENV}/bin/samtools \
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

${VENV}/lib/libz.so: ${VENV}/bin/python
	cd bin && tar -xvf zlib-1.2.8.tar.gz
	cd bin/zlib-1.2.8 && ./configure --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install
	cp -av bin/zlib-1.2.8/zlib.h ${VENV}/bin
	rm -rf bin/zlib-1.2.8

${VENV}/bin/samtools: ${VENV}/bin/python
	cd bin && tar -xvf samtools-1.3.1.tar.bz2
	cd bin/samtools-1.3.1 && ./configure --without-curses --prefix=${VIRTUAL_ENV} && $(MAKE) && $(MAKE) install

${VENV}/lib/libpng.so: ${VENV}/bin/python
	cd bin && tar -xvf libpng-1.6.24.tar.gz
	echo ${LDFLAGS}
	echo ${CFLAGS}
	echo ${LD_LIBRARY_PATH}
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
