function gen_workdir {
  export workdir=/scratch/work/bact_vcp.$jobname.$task
  mkdir $workdir
  cd $workdir
}

function finish {
  rm -rf "$workdir"
}

function install_bact_vcp {
  export PATH=$bact_vcp_path/.mc/bin:$PATH
  source activate $jobname
}

