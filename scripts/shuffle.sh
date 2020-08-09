#!/bin/bash

## shuffle fastq reads and then compress them using coreutils and pigz
## code is adapted from aaron quinlan and pierre lindenbaum on biostars
## https://www.biostars.org/p/6544/

## usage: shuffle.sh file.fastq
## usage: shuffle.sh file_1.fastq file_2.fastq

script_name=shuffle.sh

fq1=$1
fq2=$2

shuffle_paired () {
  echo Shuffling paired-end files. This may take a while ...
  echo Reading $fq1 and $fq2
  paste $fq1 $fq2 | \
  awk -v FS="\t" -v OFS="\t" '{
    head=$0
    getline seq
    getline sep
    getline qual
    print head,seq,sep,qual
  }' | \
  shuf | \
  awk -v FS="\t" -v OFS="\n" '{
    print $1,$3,$5,$7 >> ENVIRON["fq1_shuf"]
    print $2,$4,$6,$8 >> ENVIRON["fq2_shuf"]
  }' && \
  rm $fq1 $fq2
  echo Output written to $fq1_shuf and $fq2_shuf
  echo Done shuffling $fq1 and $fq2 !
}

shuffle_single () {
  echo Shuffling single-end file. This may take a while ...
  echo Reading $fq1
  awk -v FS="\t" -v OFS="\t" '{
    head=$0
    getline seq
    getline sep
    getline qual
    print head,seq,sep,qual
  }' $fq1 | \
  shuf | \
  awk -v FS="\t" -v OFS="\n" '{
    print $1,$2,$3,$4
  }' > \
  $fq_shuf && \
  rm $fq1
  echo Output written to $fq_shuf
  echo Done shuffling $fq1 !
}

compress_paired () {
  echo Compressing paired-end files. This may take a while ...
  pigz --fast -p $CPUS $fq1_shuf && mv ${fq1_shuf}.gz ${fq1}.gz
  pigz --fast -p $CPUS $fq2_shuf && mv ${fq2_shuf}.gz ${fq2}.gz
  echo Done compressing paired-end files!
  echo Reverted to original prefixes.
  echo Look for ${fq1}.gz and ${fq2}.gz
}

compress_single () {
  echo Compressing single-end file. This may take a while ...
  pigz --fast -p $CPUS $fq_shuf && mv ${fq_shuf}.gz ${fq1}.gz
  echo Done compressing single-end file!
  echo Reverted to original prefix.
  echo Look for ${fq1}.gz
}

main () {
  echo Launching $script_name
  if [ $read_format == 'paired' ]; then
    ## paired end
    export fq1_shuf=${fq1%_1.fastq}_shuf_1.fastq
    export fq2_shuf=${fq2%_2.fastq}_shuf_2.fastq
    shuffle_paired
    compress_paired
  else
    ## single end
    fq_shuf=${fq1%.fastq}_shuf.fastq
    shuffle_single
    compress_single
  fi
  echo Done with $script_name
}

main
