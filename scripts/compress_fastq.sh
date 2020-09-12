#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20000
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

target_dir=$1
CPUS=$SLURM_CPUS_PER_TASK

cd $target_dir
pigz --fast -p $CPUS *.fastq
