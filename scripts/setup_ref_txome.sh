#!/bin/bash

## -- GOAL: generate decoy-aware transcriptome
## -- first: generate herv-aware transcriptome by concatenating unannotated herv
## transcript sequences to human reference transcript sequences downloaded from
## gencode
## -- second: concatenate human reference genome (i.e. decoys) to end of
## herv-aware transcriptome
## -- third: populate decoys.txt with names of chromosomes in human reference
## genome
## -- requires config.sh to set environmental variables
## -- for use in snakemake workflow
## usage: snakemake setup_ref_txome

script_name=setup_ref_txome.sh

herv_fa=$1
base_fa=${ref_dir}assembly/gencode/gencode.v34.transcripts.fa
proj_fa=${genmod_dir}target_transcripts.fa
genome_fa=${ref_dir}assembly/encode/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
gentrome_fa=${genmod_dir}decoy_aware_transcriptome.fa
decoys=${ref_dir}assembly/encode/GRCh38_no_alt_analysis.salmon_decoys.txt

cat_hervs () {
  ## make herv-aware transcriptome
  echo Concatenating HERVs to known transcripts ...
  cat $base_fa $herv_fa > $proj_fa
  echo Output written to $proj_fa
}

cat_decoys () {
  ## make decoy-aware transcriptome
  echo Concatenating genome decoys to target transcripts ...
  cat $proj_fa $genome_fa > $gentrome_fa
  echo Output written to $gentrome_fa
}

name_decoys () {
  ## name decoys
  echo Capturing genome decoy names ...
    ## substr(s,a,b); s=string, a=start, b=end (end of string if omitted)
  awk '/^>/ {print substr($1,2)}' $genome_fa > $decoys
  echo Output written to $decoys
}

main () {
  echo Launching $script_name ...
  cat_hervs
  cat_decoys
  name_decoys
  echo Done!
}

main
