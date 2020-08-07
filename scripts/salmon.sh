#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30000
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

## -- requires config.sh to set environmental variables
## -- for use in snakemake workflow
## -- usage: snakemake salmon

## -- perform: 6h30m, 22GB for 369 gz paired samples, 8 batches of 50, 10 cpus


script_name=salmon.sh


## command-line arguments
src=$1
in_path=$2
infile_ext=$3
read_format=$4
batch_size=$5
idx_path=$6
out_path=$7
SATID=$SLURM_ARRAY_TASK_ID


startup_report () {
  ## report workflow parameters
  pwd; hostname; date
  echo Starting $script_name ...
  echo
  echo Source: $src
  echo Format: $read_format
  echo Batch size: $batch_size
  echo Input: $in_path
  echo Output: $out_path
  echo Index: $idx_path
  echo Array: $SATID
  echo
}


set_task_idx () {
  ## index code requires no "_" or "." occurs before unique identifier
  ## this is fragile and should be improved in the future

  d=$(( $batch_size * ( $SATID - 1 ) ))
  q=$(( $batch_size * $SATID ))

  if [ $read_format == 'paired' ]; then
    fs=_
  elif [ $read_format == 'single' ]; then
    fs=.
  fi

  if [ $SATID -eq 1 ]; then
    cat ${out_path}${src}.${script_name%.sh}.idx | \
      awk -v FS=$fs '{print $1}' | sort | uniq | sed -e "${q}q" > \
      ${out_path}${src}_${script_name%.sh}_idx_${SATID}.tmp
  elif [ $SATID -gt 1 ]; then
    cat ${out_path}${src}.${script_name%.sh}.idx | \
      awk -v FS=$fs '{print $1}' | sort | uniq | sed -e "1,${d}d;${q}q" > \
      ${out_path}${src}_${script_name%.sh}_idx_${SATID}.tmp
  fi
}


reset_path () {
  ## purge non-default paths
  ## add paths to required binaries
  echo Resetting \$PATH ...
  module list
  echo Currently loaded modules will be purged
  module purge
  echo Adding gcc/7.1.0 and salmon/1.2.1 to \$PATH ...
  module load gcc/7.1.0 salmon/1.2.1  ## this adds salmon bin to PATH
}


shutdown_report () {
  date
  echo end $script_name
  echo seff $SLURM_JOBID
  echo sacct -o reqmem,maxrss,elapsed,alloccpus,nodelist -j $SLURM_JOBID
}


run_module () {
  ## initiation workflow loop
  for elem in $(<${out_path}${src}_${script_name%.sh}_idx_${SATID}.tmp); do  ## elem includes path to input
    sn=`echo $elem | awk -v FS='/' '{print $NF}'`
    echo Starting quantification for sample $sn

    ## if src == barcuva, reads are paired, length is 51 or 101
    ## if src == tcga, sra some reads are paired and others are single,
    ##    length is variable

    if [ $src == 'barcuva' ]; then
      salmon quant \
        --gcBias \
        -p $SLURM_CPUS_PER_TASK \
        -i $idx_path \
        -l A \
        -1 <(gunzip -c ${elem}_R1_val_1${infile_ext}) \
        -2 <(gunzip -c ${elem}_R2_val_2${infile_ext}) \
        --validateMappings \
        -o ${out_path}${sn}
    else
      if [ $read_format == 'paired' ]; then
        salmon quant \
          --gcBias \
          -p $SLURM_CPUS_PER_TASK \
          -i $idx_path \
          -l A \
          -1 <(gunzip -c ${elem}_1${infile_ext}) \
          -2 <(gunzip -c ${elem}_2${infile_ext}) \
          --validateMappings \
          -o ${out_path}${sn}
      elif [ $read_format == 'single' ]; then
        salmon quant \
          --gcBias \
          -p $SLURM_CPUS_PER_TASK \
          -i $idx_path \
          -l A \
          -r <(gunzip -c ${elem}${infile_ext}) \
          --validateMappings \
          -o ${out_path}${sn}
      else
        echo Encountered a problem with $sn
        exit
      fi
    fi
    echo Done with quantification for sample $sn
  done
  rm ${out_path}${src}_${script_name%.sh}_idx_${SATID}.tmp
}


main () {
  startup_report
  set_task_idx
  reset_path
  run_module
  shutdown_report
}


main
