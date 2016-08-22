function gen_workdir {
  export workdir=/scratch/work/bact_vcp.$jobname.$task
  mkdir $workdir
  cd $workdir
}

function finish {
  rm -rf "$workdir"
}

function install_bact_vcp {
  cp -a $bact_vcp_path $workdir
  export PATH=$workdir/.mc/bin:$PATH
  source activate venv
}



