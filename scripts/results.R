## R code to show results of differential gene expression analysis

## load field effect results and full herv results into dataframes
## display pearson correlations for cohort A, NATvHLT, TUMvHLT, TUMvNAT
## intersection will capture protein coding genes only

## count total herv transcripts, select herv transcripts, select herv genes,
## and moderately expressed herv genes and display numbers in bar chart
## requires attig gtf, herv_ids, and dgeobj$dds

## display differentially expressed herv genes on background of protein coding
## genes with a volcano plot (L2FC vs -log10 pval)
## requires dgeobj

## load dgeobj into list and herv-specific results into dataframe
## create boxplots showing top herv hits

## try to highlight samples with MSS

## try showing alignments of top hit along with sequence and predictied protein
## try showing possible HLA binding sites in predicted protein


BiocManager::install()
library(readr)
library(DESeq2)
library(ashr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(bacon)


args <- commandArgs(trailingOnly=TRUE)
cohort <- args[2]
cat(paste("Cohort:", cohort, "\n"))


source("res_utilities.R")
source("res_themes.R")
source("res_corr.R")
source("res_hervsum.R")
source("res_volcano.R")
source("res_topherv.R")


correlation_analysis <- function() {
    tab_herv <- load_res_tables("herv")
    tab_field <- load_res_tables("field")
    l <- calc_cor(tab_herv, tab_field)
    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-field-cor_A.pdf")
    write_cor(l$fig, target)
}


herv_element_summary <- function() {
    gtf <- paste0(Sys.getenv("attig_dir"), "attig_flat_transcriptome.gtf.gz")
    ids <- paste0(Sys.getenv("genmod_dir"), "herv_ids.tsv")
    hrvnums <- count_herv_ids(gtf, ids)
    target <- paste0(Sys.getenv("rdata_dir"), "dge_", cohort, ".Rda")
    dgeobj <- load_dgeobj(target)
    hrvnums[["Expressed HERV Genes"]] <- sum(grepl("HERV", rownames(dgeobj$dds)))
    hrvnames <- rownames(dgeobj$dds)[grepl("HERV", rownames(dgeobj$dds))]
    hrvcats <- count_herv_cats(hrvnames)
    sumlist <- setNames(list(hrvnums, hrvcats), c("HERV Genes", "HERV Classes"))
    pl <- plot_hrvsum(sumlist)
    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-element-sum_A.pdf")
    write_hrvsum(pl, target)
    return(dgeobj)
}


volcano_analysis <- function(dgeobj) {
    pl <- assemble_volcano(dgeobj)
    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-volcano_", cohort, ".pdf")
    write_volcano(pl, target)
}


topherv_analysis <- function(dgeobj) {
    tophervs <- select_tophervs()
    dfs <- format_counts(tophervs, dgeobj)
    pl <- plot_counts(dfs)
    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-counts_", cohort, ".pdf")
    write_counts(pl, target)
}


main <- function() {
    if (cohort == "A") {
        correlation_analysis()
        dgeobj <- herv_element_summary()
    } else {
        target <- paste0(Sys.getenv("rdata_dir"), "dge_", cohort, ".Rda")
        dgeobj <- load_dgeobj(target)
    }
    volcano_analysis(dgeobj)
    if (cohort %in% c("A", "C")) {
        topherv_analysis(dgeobj)
    }
}


main()




# target <- "/scratch/chd5n/test-fig.pdf"
# ggsave(target, p, device="pdf", width=15, height=4, units="in")

tab_herv <- load_res_tables("herv")
target_dir <- Sys.getenv("plot_dir")
check_dir(target_dir)
target <- paste0(target_dir, "herv-qq_", cohort, ".pdf")
simple_qq(tab_herv, target)
