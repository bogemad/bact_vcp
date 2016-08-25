#!/bin/bash
#PBS -N 
#PBS -l ncpus=1
#PBS -l mem=4gb
#PBS -l walltime=4:00:00
#PBS -q smallq
#PBS -o 
#PBS -e 

jobname=
task=
bact_vcp_path=
source $bact_vcp_path/scripts/hpc.job_lib.sh
gen_workdir
trap finish EXIT
install_bact_vcp

