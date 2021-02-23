## find colorectal RNA-seq in ORIEN

set_vars ()
{
    data_root=/project/orien/data/
    target=${data_root}Linkage_and_QC_Files/UVA_Clinical_Specimen_Linkage_Data_Full_20200901.csv
    herv_dir=/scratch/chd5n/herv/
}

define_samples ()
{
    awk -v FS="," '
        NR == 1 || $8 ~ /GI - Colorectal Cancer/ {print $0}
    ' ${target} > ${herv_dir}orien_crc_samples.csv

    find ${data_root}Avatar_MolecularData_hg38/ -type f -name "*.fastq.gz" \
    | grep 'RNA' > ${herv_dir}orien_rna_fastq.txt

    awk -v FS="," '
        NR == FNR {gsub(/"/, "", $6); a[$6]; next}
        {split($0, l, "/"); split(l[9], l2, "_"); if (l2[1] in a) print $0}
    ' ${herv_dir}orien_crc_samples.csv ${herv_dir}orien_rna_fastq.txt \
    > ${herv_dir}orien_crc_fastq.txt

    printf "%s\n" "specimen data written to ${herv_dir}orien_crc_samples.csv"
    printf "%s\n" "fastq inventory written to ${herv_dir}orien_rna_fastq.txt"
    printf "%s\n" "crc fastq paths written to ${herv_dir}orien_crc_fastq.txt"
}

main ()
{
    set_vars
    define_samples
}

main
