
VENV = .venv
export VIRTUAL_ENV := $(abspath ${VENV})
export PATH := ${VIRTUAL_ENV}/bin:${PATH}


all: bin/smalt-0.7.6/src/smalt \
${VENV}/bin/python \
${VENV}/bin/cutadapt \

clean: 
	rm -rf bin/smalt-0.7.6 ${VENV} intermediate_files

.PHONY: all clean
.SECONDARY:

bin/smalt-0.7.6/src/smalt:
	cd bin && tar -xvf smalt-0.7.6-static.tar.gz
	cd bin/smalt-0.7.6 && ./configure && $(MAKE)


${VENV}/bin/python:
	python2 support/virtualenv.py ${VENV}

${KEYS} ${RESULTS} ${INT_FILES}:
	mkdir $@

${VENV}/bin/cutadapt: ${VENV}/bin/python
	pip install --upgrade -r scripts/requirements.pip
