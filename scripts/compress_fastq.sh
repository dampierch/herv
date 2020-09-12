#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=01000
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

## 8MB, 1h45m for 59 fastq samples, paired, single mix

target_dir=$1
CPUS=$SLURM_CPUS_PER_TASK

cd $target_dir
pigz --fast -p $CPUS *.fastq

echo seff $SLURM_JOBID
echo sacct -o reqmem,maxrss,elapsed,alloccpus,nodelist -j $SLURM_JOBID
