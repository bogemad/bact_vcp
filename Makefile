all: make_smalt make_samtools ${VENV} python-reqs

clean: 
	rm -rf smalt samtools

.PHONY: all clean
.SECONDARY:

make_smalt: smalt/configure
	cd smalt && ./configure && $(MAKE)

make_samtools: samtools/configure
	cd samtools && ./configure && $(MAKE)

VENV = .venv
export VIRTUAL_ENV := $(abspath ${VENV})
export PATH := ${VIRTUAL_ENV}/bin:${PATH}

${VENV}:
	python3 -m venv $@

python-reqs: requirements.pip | ${VENV}
	pip install --upgrade -r requirements.pip

