## R code to perform differential gene expression analysis

## load txi into list and pheno file into dataframe
## create DESeq2 data set with txi and dataframe
## pre-filter genes prior to normalization
## get size-factors and normalized counts
## prepare input for sva and estimate number of latent factors
## estimate latent factors themselves
## adjust counts for batch
## filter for protein coding and herv genes
## test for differential expression, extract and display results


BiocManager::install()
library(readr)
library(DESeq2)
library(sva)
library(limma)
library(ensembldb)
library(EnsDb.Hsapiens.v86)


args <- commandArgs(trailingOnly=TRUE)
cohort <- args[2]
cat(paste("Cohort:", cohort, "\n"))


load_inputs <- function(target1, target2, target3) {
    ## load txi into list and pheno file into dataframe
    cat("Loading inputs\n")
    nm <- base::load(target1, verbose=TRUE)
    df <- data.frame(
        readr::read_tsv(target2)
    )
    df$phenotype <- factor(df$phenotype, levels=c("HLT", "NAT", "CRC"))
    hv <- data.frame(
        readr::read_tsv(target3)
    )
    v <- as.vector(hv$herv_id)
    return(setNames(list(txi, df, v), c("txi", "colData", "hervs")))
}


create_dds <- function(l) {
    ## create DESeq2 data set with txi and dataframe
    ## pre-filter genes prior to normalization
    ## get size-factors and normalized counts
    cat("Creating pre-filtered dds\n")
    dds <- DESeqDataSetFromTximport(l$txi, l$colData, ~ phenotype)
    count_lim <- 0
    sample_lim <- 1/2 * ncol(dds)
    keep <- rowSums( counts(dds) > count_lim ) >= sample_lim
    dds <- dds[keep, ]
    dds <- DESeq(dds)
    return(dds)
}


prep_sva <- function(dds) {
    ## prepare input for sva and estimate number of latent factors
    cat("Preparing data for SVA\n")
    dat <- counts(dds, normalized=TRUE)  ## extract the normalized counts
    mod <- model.matrix(~phenotype, data=colData(dds))  ## set the full model
    mod0 <- model.matrix(~1, data=colData(dds))  ## set the null model
    nsv <- num.sv(dat, mod, method=c("be"), B=20, seed=1)  ## Buja Eyuboglu method
    cat(nsv, "significant latent factors\n")
    svaobj <- setNames(
        list(dat, mod, mod0, nsv),
        c("counts", "full", "null", "nsv")
    )
    return(svaobj)
}


run_sva <-function(l, dds) {
    ## estimate latent factors themselves
    ## append latent factors to dds
    cat("Estimating surrogate variables\n")
    svs <- svaseq(l$counts, l$full, l$null, n.sv=l$nsv)
    cat("\nAdding surrogate variables to dds\n")
    for (i in seq_len(svs$n.sv)) {
        newvar <- paste0("sv", i)
        colData(dds)[, newvar] <- svs$sv[, i]
    }
    nvidx <- (ncol(colData(dds)) - i + 1):ncol(colData(dds))
    newvars <- colnames(colData(dds))[nvidx]
    d <- formula(
      paste0("~", paste(paste(newvars, collapse="+"), "phenotype", sep="+"))
    )
    design(dds) <- d
    return(dds)
}

batch_adjustment <- function(dds) {
    ## adjust counts for batch
    cat("Adjusting counts for batch\n")
    vsd <- vst(dds, blind=FALSE)
    counts <- assay(vsd)
    design <- model.matrix(~ phenotype, data=colData(vsd))
    idx <- base::grepl("sv\\d", colnames(colData(dds)))
    covs <- colnames(colData(dds))[idx]
    covariates <- colData(vsd)[, covs]
    correctedY <- limma::removeBatchEffect(
      counts, covariates=covariates, design=design
    )
    correctedY <- t(correctedY)
    correctedY[correctedY < 0] <- 0
    adj_counts <- t(correctedY)
    return(adj_counts)
}

filter_genes <- function(dds, adj, hervs) {
    ## filter for protein coding and moderately expressed herv genes
    cat("Filtering genes for test set\n")
    df <- data.frame(
        ensembldb::genes(EnsDb.Hsapiens.v86)
    )
    idx <- df$gene_biotype == "protein_coding"
    pcg <- as.vector(df[idx, "gene_id"])
    idx <- rownames(adj) %in% pcg
    adjpcg <- adj[idx, ]
    idx <- rownames(adj) %in% hervs
    adjhrv <- adj[idx, ]
    targets <- pcg
    for (e in levels(colData(dds)$phenotype)) {
        mdn <- median(adjpcg[, colData(dds)$phenotype == e])
        std <- sd(adjpcg[, colData(dds)$phenotype == e])
        cutoff <- mdn - std
        keep <- apply(adjhrv[, colData(dds)$phenotype == e], 1, median) > cutoff
        targets <- c(targets, rownames(adjhrv)[keep])
    }
    test_set <- base::unique(targets)
    idx <- rownames(dds) %in% test_set
    dds <- dds[idx, ]
    cat("Returning", nrow(dds), "genes in test set\n")
    return(dds)
}

test_dge <- function(dds) {
    ## test for differential expression, extract and display results
    cat("Testing for differential expression\n")
    dds <- DESeq(dds)
    cons <- list()
    m <- combn(levels(colData(dds)$phenotype), 2)
    for (i in seq_len(ncol(m))) {
        cons[[i]] <- c("phenotype", rev(m[, i]))
        names(cons) <- c(
            names(cons)[seq_len(length(cons) - 1)], paste(rev(m[, i]), collapse="v")
        )
    }
    res <- list()
    for (i in seq_len(length(cons))) {
        res[[i]] <- results(dds, contrast=cons[[i]], alpha=0.05)  ## default alpha is 0.1
    }
    names(res) <- names(cons)
    return(setNames(list(dds, res), c("final", "results")))
}

write_dge <- function(dgeobj, target) {
  save(dgeobj, file=target)
  cat("dgeobj written to", target, "\n")
}

main <- function() {
    target1 <- paste0(Sys.getenv("rdata_dir"), "txi_", cohort, ".Rda")
    target2 <- paste0(Sys.getenv("ann_dir"), "pheno_", cohort, ".tsv")
    target3 <- paste0(Sys.getenv("genmod_dir"), "herv_ids.tsv")
    inputs <- load_inputs(target1, target2, target3)
    dds <- create_dds(inputs)
    svaobj <- prep_sva(dds)
    ddssva <- run_sva(svaobj, dds)
    vsdadj <- batch_adjustment(ddssva)
    ddssva <- filter_genes(ddssva, vsdadj, inputs$hervs)
    dgeobj <- test_dge(ddssva)
    target <- paste0(Sys.getenv("rdata_dir"), "dge_", cohort, ".Rda")
    write_dge(dgeobj, target)

    print(lapply(dgeobj$results, head))
    print(lapply(dgeobj$results, summary))
}

main()

# x <- dgeobj$results$CRCvNAT[order(dgeobj$results$CRCvNAT[,"padj"], decreasing=FALSE),]


# r <- res[[1]][order(res[[1]]$padj), ]
# r <- res[[2]][order(res[[2]]$padj), ]
# r <- res[[3]][order(res[[3]]$padj), ]
# r$rank <- seq(1, nrow(r))
#
# ## save gene lists
# setwd("/scratch/chd5n/")
# file <- paste0("natvhlt.tsv")
# file <- paste0("crcvhlt.tsv")
# file <- paste0("crcvnat.tsv")
# write.table(r, file=file, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
