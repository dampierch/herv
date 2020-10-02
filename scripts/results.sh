#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=02:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

## -- requires config.sh to set environmental variables
## -- for use in snakemake workflow
## -- usage: snakemake setup_dge

## -- perform: 0.5h, 10GB for cohort A


script_name=results.sh


## command-line arguments
cohort=$1


startup_report () {
  ## report workflow parameters
  pwd; hostname; date
  echo Starting $script_name ...
  echo
  echo Cohort: $cohort
  echo
}


reset_path () {
  ## purge non-default paths
  ## add paths to required binaries
  echo Resetting \$PATH ...
  module list
  echo Currently loaded modules will be purged
  module purge
  echo Adding gcc/7.1.0 openmpi/3.1.4 R/4.0.0 to \$PATH ...
  module load gcc/7.1.0 openmpi/3.1.4 R/4.0.0  ## this adds R bin to PATH
}


shutdown_report () {
  date
  echo end $script_name
  echo seff $SLURM_JOBID
  echo sacct -o reqmem,maxrss,elapsed,alloccpus,nodelist -j $SLURM_JOBID
}


run_module () {
  echo Running R script ...
  Rscript ${script_name%.sh}.R --args $cohort
}


main () {
  startup_report
  reset_path
  run_module
  shutdown_report
}


main
