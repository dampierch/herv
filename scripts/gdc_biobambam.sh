#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50000
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

## -- requires config.sh to set environmental variables
## -- for use in snakemake workflow
## -- usage: snakemake convert_tcga

## -- perform: 1h, 41MB for 220 bam single, 8 batches of 30, 1 cpu w/o shuffle
## -- perform: 2h20m, 8GB for 220 bam single, 8 batches of 30, 4 cpu w/ shuffle


script_name=gdc_biobambam.sh


## command-line arguments
in_path=$1
infile_ext=$2
export read_format=$3  ## shuffle needs this variable
batch_size=$4
out_path=$5
SATID=$SLURM_ARRAY_TASK_ID
export CPUS=$SLURM_CPUS_PER_TASK  ## shuffle needs this variable for pigz


startup_report () {
  ## report workflow parameters
  pwd; hostname; date
  echo Starting $script_name ...
  echo
  echo Format: $read_format
  echo Batch size: $batch_size
  echo Input: $in_path
  echo Output: $out_path
  echo Array: $SATID
  echo CPUs: $CPUS
  echo
}


set_task_idx () {
  ## single input per sample obviates need for selection of unique samples

  d=$(( $batch_size * ( $SATID - 1 ) ))
  q=$(( $batch_size * $SATID ))

  if [ $SATID -eq 1 ]; then
    cat ${out_path}${script_name%.sh}.idx | \
      sed -e "${q}q" > ${out_path}${script_name%.sh}_idx_${SATID}.tmp
  elif [ $SATID -gt 1 ]; then
    cat ${out_path}${script_name%.sh}.idx | \
      sed -e "1,${d}d;${q}q" > ${out_path}${script_name%.sh}_idx_${SATID}.tmp
  fi
}


reset_path () {
  ## purge non-default paths
  ## add paths to required binaries
  echo Resetting \$PATH ...
  module list
  echo Currently loaded modules will be purged
  module purge
  echo Adding biobambam2/2.0.87 to \$PATH ...
  PATH=$PATH:~/apps/biobambam2/2.0.87-release-20180301132713/x86_64-etch-linux-gnu/bin/  ## this adds biobambam bin to PATH
}


shutdown_report () {
  date
  echo end $script_name
  echo seff $SLURM_JOBID
  echo sacct -o reqmem,maxrss,elapsed,alloccpus,nodelist -j $SLURM_JOBID
}


run_module () {
  ## initiation workflow loop
  for elem in $(<${out_path}${script_name%.sh}_idx_${SATID}.tmp); do  ## elem includes path to input
    sn=`echo $elem | awk -v FS="::" '{print $1}'`
    fn=`echo $elem | awk -v FS="::" '{print $2}'`
    echo Starting conversion for sample $sn
    echo '--' File name: $fn

    ## tcga, so some reads are paired and others are single, length is variable

    if [ $read_format == 'paired' ]; then
      bamtofastq \
        gz=0 \
        filename=${fn} \
        F=${out_path}${sn}_1.fastq \
        F2=${out_path}${sn}_2.fastq
      ${command_dir}shuffle.sh ${out_path}${sn}_1.fastq ${out_path}${sn}_2.fastq
    elif [ $read_format == 'single' ]; then
      bamtofastq \
        gz=0 \
        filename=${fn} \
        S=${out_path}${sn}.fastq
      ${command_dir}shuffle.sh ${out_path}${sn}.fastq
    else
      echo Encountered a problem with $sn
      exit
    fi
    echo Done with conversion for sample $sn
  done
  rm ${out_path}${script_name%.sh}_idx_${SATID}.tmp
}


main () {
  startup_report
  set_task_idx
  reset_path
  run_module
  shutdown_report
}


main
