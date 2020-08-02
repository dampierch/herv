#!/bin/bash

## config file for project
## to get these variables into parent shell, use source command

script_name=config.sh


proj_name=herv


export_home_dirs () {
  export proj_home=${HOME}/projects/${proj_name}/
  export command_dir=${proj_home}scripts/
  export log_dir=${proj_home}logs/
}

export_data_dirs () {
  export data_root=/scratch/${USER}/
  export cphg_root=/project/CPHG/CASEY/
  export work_dir=${data_root}${proj_name}/
  export ref_dir=${data_root}reference-genome/
  export download_dir=${data_root}downloads/
}

export_proj_top_dirs () {
  export ann_dir=${work_dir}annotations/
  export fastq_dir=${work_dir}fastq/
  export fqc_dir=${work_dir}fqc/
  export mqc_dir=${work_dir}mqc/
  export star_dir=${work_dir}star/
  export rsem_dir=${work_dir}rsem/
  export salmon_dir=${work_dir}salmon/
  export rdata_dir=${work_dir}rdata/
  export res_dir=${work_dir}results/
}

export_proj_sub_dirs () {
  export plot_dir=${res_dir}plots/
  export table_dir=${res_dir}tables/
}

export_proj_special_dirs () {
  export attig_dir=${work_dir}attig/
  export genmod_dir=${work_dir}genmod/
  export barcuva_fqdir=${cphg_root}chris-misc/fastq/
  export salidx_dir=${salmon_dir}gentrome_k25/
  export salqnt_dir=${salmon_dir}quant/
}

check_directories () {
  home_list="proj_home command_dir log_dir"
  data_list="data_root cphg_root work_dir ref_dir download_dir"
  proj_list_top="ann_dir fastq_dir fqc_dir mqc_dir star_dir rsem_dir salmon_dir rdata_dir res_dir"
  proj_list_sub="plot_dir table_dir"
  special_list="attig_dir genmod_dir barcuva_fqdir salidx_dir salqnt_dir"

  for e in $home_list $data_list $proj_list_top $proj_list_sub $special_list; do
    if [[ -e ${!e} ]]; then
      echo ${!e} present
    else
      echo ${!e} absent
    fi
  done
}

main () {
  echo Launching $script_name for project $proj_name ...
  echo Exporting variables to shell environment ...
  export_home_dirs
  export_data_dirs
  export_proj_top_dirs
  export_proj_sub_dirs
  export_proj_special_dirs
  echo Checking for presence/absence of directories ...
  check_directories
  echo Done!
  echo Tip: Use source command if the variables you were expecting cannot be seen in parent shell.
}

main
