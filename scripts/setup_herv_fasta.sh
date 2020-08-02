#!/bin/bash

## -- input reference genome in FASTA and HERV transcriptome in BED12, extract
## HERV transcript sequences in FASTA
## -- input must be unzipped
## -- requires config.sh to set environmental variables
## -- for use in snakemake workflow
## usage: snakemake setup_herv_fasta

script_name=setup_herv_fasta.sh

bed=$1
ref=${ref_dir}assembly/encode/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
fasta=${genmod_dir}herv_transcripts.fa

reset_paths () {
  ## purge non-default paths
  ## add paths to required binaries
  echo Resetting \$PATH ...
  module list
  echo Currently loaded modules will be purged
  module purge
  echo Adding gcc/7.1.0 and bedtools/2.26.0 to \$PATH ...
  module load gcc/7.1.0 bedtools/2.26.0  ## this adds bedtools to PATH
}

extract_seq () {
  ## extracts target sequences defined in BED from reference FASTA
  echo Reading reference FASTA from $ref
  echo Reading target BED from $bed
  bedtools getfasta \
    -split \
    -name \
    -fi $ref \
    -bed $bed | \
    fold -w 60 > \
    $fasta
  echo Output FASTA written to $fasta
}

main() {
  echo Launching $script_name
  reset_paths
  echo bedtools getfasta is working ...
  extract_seq
  echo Done!
}

main
