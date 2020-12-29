# Quick Check on Protein Coding Potential

## Xp22 Locus

```
cd /scratch/chd5n/herv/genmod/
grep 'LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320' herv_ids.tsv > Xp22_ids.tsv
```

The Xp22 locus harbors the following sequences:

transbcbde997974f0776::chrX:4539292-4543738     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
*transec049363fddb03a2::chrX:4540474-4541765     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
trans6dd9b80376bea865::chrX:4540514-4546155     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
transd71a404b5b4e5332::chrX:4540659-4546323     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
transf1650e45563ccea6::chrX:4540890-4545238     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
trans950d281b6b8a110b::chrX:4540919-4542591     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
trans04b631de9ba0dcaa::chrX:4541445-4542986     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
*trans49fae8cab1376dbb::chrX:4541537-4545907     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
transc1b3baa858d598f6::chrX:4543716-4546070     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
trans1c69e5d603801f57::chrX:4544781-4545355     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
transa81aa55c07885ccb::chrX:4544801-4545462     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
transbc2cb9416cf409ac::chrX:4544920-4545528     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320
trans786c3b072fa4b5ce::chrX:4545046-4546135     LTR/ERV1|HERVH-int~LTR7Y|X|4540474|4546320

After checking counts for each of the 13 transcripts in two random TCGA tumor quant.sf files, we find the two transcripts with the highest expression are:

1. transec049363fddb03a2::chrX:4540474-4541765
2. trans49fae8cab1376dbb::chrX:4541537-4545907

In order to check protein homology for the predicted proteins that could be translated from those transcripts, we first have to extract their sequences. To do that, we will use BEDTools, which requires a BED as input.

```
grep 'trans49fae8cab1376dbb' herv_transcriptome.bed > Xp22.bed
grep 'transec049363fddb03a2' herv_transcriptome.bed >> Xp22.bed
```

This is what the BED looks like for the two highly expressed transcripts from the Xp22 locus:

```
chrX	4541537	4545907	trans49fae8cab1376dbb	.	.	.	.	.	2	3485,691	0,3679
chrX	4540474	4541765	transec049363fddb03a2	.	.	.	.	.	1	1291	0
```

We make a FASTA for them with BEDTools.

```
module load gcc/7.1.0 bedtools/2.26.0

bed=Xp22.bed
ref=/scratch/chd5n/reference-genome/assembly/encode/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
fasta=Xp22.fa

bedtools getfasta \
    -split \
    -name \
    -fi $ref \
    -bed $bed | \
    fold -w 60 > \
    $fasta
```

We open the FASTA, manually copy and paste one of the sequences into BLAST (via website), and find some reasonable homologous proteins and exactly the right DNA alignments. Doing this systematically will require protein expertise and some BLAST scripts. That is for the next project.
