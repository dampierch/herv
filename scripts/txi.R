## R code to make txi objects with tximport

BiocManager::install()
library(tximport)
library(readr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

args <- commandArgs(trailingOnly=TRUE)
cohort <- args[2]
cat(paste("Cohort:", cohort, "\n"))

setup_tx_map <- function(target) {
  ## create tx2gene dataframe for tximport (2 columns: txID, geneID)
  hvdf <- data.frame(
    readr::read_tsv(target)
  )
  colnames(hvdf) <- c("tx_id", "gene_id")
  txdf <- data.frame(
    ensembldb::transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_id"))
  )
  txdf <- subset(txdf, select=c(tx_id, gene_id))
  df <- rbind(txdf, hvdf)
  return(df)
}

write_txi <- function(txi, target) {
  save(txi, file=target)
  cat("txi written to", target, "\n")
}

main <- function() {
  target <- paste0(Sys.getenv("genmod_dir"), "herv_ids.tsv")
  txmap <- setup_tx_map(target)
  target <- paste0(Sys.getenv("ann_dir"), "pheno_", cohort, ".tsv")
  p <- readr::read_tsv(target)
  targets <- p$path
  txi <- tximport(targets, type="salmon", tx2gene=txmap, ignoreTxVersion=TRUE)
  target <- paste0(Sys.getenv("rdata_dir"), "txi_", cohort, ".Rda")
  write_txi(txi, target)
}

main()
