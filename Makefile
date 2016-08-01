all: make_smalt make_ngsutils ${VENV} python-reqs 

clean: 


.PHONY: all clean
.SECONDARY:

make_smalt:
	cd bin/smalt && ./configure && $(MAKE)

make_ngsutils:
	cd bin/ngsutils && $(MAKE)

VENV = .venv
export VIRTUAL_ENV := $(abspath ${VENV})
export PATH := ${VIRTUAL_ENV}/bin:${PATH}

${VENV}: 
	python2 -m virtualenv $@

python-reqs: scripts/requirements.pip | ${VENV}
	pip install --upgrade -r scripts/requirements.pip

