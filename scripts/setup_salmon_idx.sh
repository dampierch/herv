#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=30000
#SBATCH --time=1:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH	--mail-type=END
#SBATCH	--mail-user=chd5n@virginia.edu


## -- this script assumes a reference genome and transcriptome are already
## downloaded and concatenated into a reference 'gentrome.fa' file and that
## gentrome decoy names are listed in decoys.txt file
## -- generates salmon index for quantification with selective alignment
## -- requires config.sh to set environmental variables
## -- for use in snakemake workflow
## usage: snakemake setup_salmon_idx

## -- takes about 21 GB, 25 min for human on 12 cores

script_name=setup_salmon_idx.sh

## header
pwd; hostname; date

## set variables
gentrome=$1
decoys=$2
kmer=$3  ## 31 is appropriate default for 75bp and longer reads
cpus_per_task=${SLURM_CPUS_PER_TASK}

prep_directory () {
  ## prepare salmon index directory structure
  if [[ -e ${salmon_dir} ]]; then
    :
  else
    mkdir ${salmon_dir}
  fi

  if [[ -e ${salmon_dir}gentrome_k${kmer} ]]; then
    :
  else
    mkdir ${salmon_dir}gentrome_k${kmer}
  fi
}

reset_paths () {
  ## purge non-default paths
  ## add paths to required binaries
  echo Resetting \$PATH ...
  module list
  echo Currently loaded modules will be purged
  module purge
  echo Adding gcc/7.1.0 and salmon/1.2.1 to \$PATH ...
  module load gcc/7.1.0 salmon/1.2.1  ## this adds salmon bin to PATH
}

make_index () {
  ## generate salmon index for given k-mer length
  echo Reading gentrome from $gentrome ...
  echo Reading decoys from $decoys ...
  salmon index \
    -p $cpus_per_task \
    -t $gentrome \
    -d $decoys \
    -i ${salmon_dir}gentrome_k${kmer} \
    -k $kmer \
    --gencode
  echo Output written to ${salmon_dir}gentrome_k${kmer}
}

main () {
  echo Launching $script_name ...
  echo Preparing directory structure for output ...
  prep_directory
  reset_paths
  make_index
  echo Done with ${script_name}!
  date
  echo seff $SLURM_JOBID
  echo sacct -o reqmem,maxrss,elapsed,alloccpus,nodelist -j $SLURM_JOBID
}

main
