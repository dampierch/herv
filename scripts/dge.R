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
    if (cohort != "B") {
        df$phenotype <- factor(df$phenotype, levels=c("HLT", "NAT", "CRC"))
    } else {
        df$phenotype <- factor(df$phenotype, levels=c("NAT", "CRC"))
        df$subject_id <- factor(df$subject_id)
    }
    hv <- data.frame(
        readr::read_tsv(target3)
    )
    v <- as.vector(hv$herv_id)
    return(setNames(list(txi, df, v), c("txi", "colData", "hervs")))
}


create_dds <- function(l) {
    ## create DESeq2 data set with txi and dataframe
    ## pre-filter genes prior to normalization
    ## get normalizationFactors for normalized counts
    cat("Creating pre-filtered dds\n")
    if (cohort != "B") {
        dds <- DESeqDataSetFromTximport(l$txi, l$colData, ~ phenotype)
    } else {
        dds <- DESeqDataSetFromTximport(l$txi, l$colData, ~ subject_id + phenotype)
    }
    count_lim <- 0
    sample_lim <- 1/2 * ncol(dds)
    keep <- rowSums( counts(dds) > count_lim ) >= sample_lim
    dds <- dds[keep, ]
    dds <- estimateSizeFactors(dds)
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


count_adjustment <- function(dds) {
    ## adjust counts with variance stabilizing transformation for all cohorts
    ## adjust counts for batch effect for cohorts A and C
    cat("Adjusting counts\n")
    vsd <- vst(dds, blind=FALSE)
    counts <- assay(vsd)
    if (cohort != "B") {
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
    } else {
        adj_counts <- counts
    }
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
    idx <- rownames(adj) %in% test_set
    adj <- adj[idx, ]
    cat("Returning", nrow(dds), "genes in test set\n")
    return(setNames(list(dds, adj), c("dds", "adj")))
}


map_ids <- function(res) {
    ## mapping of ensembl gene ids to gene symbols and entrez ids
    tab <- subset(
        data.frame(
            ensembldb::genes(
                EnsDb.Hsapiens.v86,
                columns=c("gene_id", "gene_name", "entrezid")
            )
        ),
        select=c("gene_id", "gene_name", "entrezid")
    )
    x <- res[, "gene_id"]
    idx <- match(x, tab$gene_id)
    res[, "gene_name"] <- tab[idx, "gene_name"]
    res[, "entrez_id"] <- tab[idx, "entrezid"]
    res <- res[, c(8:10,1:7)]
    return(res)
}


test_dge <- function(dds, adj, fc=2, alpha=0.05) {
    ## test for differential expression, extract and sort results
    ## adj counts added to dgeobj for downstream visualization
    cat("Testing for differential expression\n")
    dds <- DESeq(dds)
    cons <- list()
    m <- combn(levels(colData(dds)$phenotype), 2)
    for (i in seq_len(ncol(m))) {
        cons[[i]] <- c("phenotype", rev(m[, i]))
        names(cons) <- c(
            names(cons)[seq_len(length(cons) - 1)],
            paste(rev(m[, i]), collapse="v")
        )
    }
    res <- list()
    for (i in seq_len(length(cons))) {
        res[[i]] <- results(
            dds,
            lfcThreshold=log2(fc),
            altHypothesis="greaterAbs",
            contrast=cons[[i]],
            alpha=alpha
        )
        res[[i]] <- res[[i]][order(res[[i]][, "padj"], decreasing=FALSE), ]
        res[[i]][, "rank"] <- seq(1, nrow(res[[i]]))
        res[[i]][, "gene_id"] <- rownames(res[[i]])
        res[[i]] <- map_ids(res[[i]])
    }
    names(res) <- names(cons)
    return(setNames(list(dds, adj, res), c("dds", "adj", "res")))
}


check_dir <- function(dir) {
    ## check if intended output directory exists and create if it does not
    d <- "/"
    for (e in unlist(strsplit(dir, "/"))[2:length(unlist(strsplit(dir, "/")))]) {
        d <- paste0(d, e, "/")
        if (file.exists(d)) {
            cat(d, "exists\n")
        } else {
            cat(d, "does not exist; creating now\n")
            dir.create(d)
        }
    }
}


write_dds <- function(ddsobj, target) {
    ## save prelim dds and adjusted counts
    cat("Writing ddsobj to file\n")
    save(ddsobj, file=target)
    cat("ddsobj written to", target, "\n")
}


write_dge <- function(dgeobj, target) {
    cat("Writing dgeobj to file\n")
    save(dgeobj, file=target)
    cat("dgeobj written to", target, "\n")
}


write_gene_list <- function(res, herv=FALSE) {
    cat("Writing results table to file\n")
    for (i in seq_len(length(res))) {
        con <- names(res)[i]
        if (herv) {
            target <- paste0(
                Sys.getenv("table_dir"), "res_", cohort, "_", con, "_hrv.tsv"
            )
        } else {
            target <- paste0(
                Sys.getenv("table_dir"), "res_", cohort, "_", con, "_all.tsv"
            )
        }
        write.table(res[[i]], file=target, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=TRUE)
        cat("gene list written to", target, "\n")
    }
}


extract_hervs <- function(res) {
    ## extract herv-specific results table
    cat("Extracting HERV-specific results\n")
    for (i in seq_len(length(res))) {
        keep <- grepl("HERV", res[[i]][, "gene_id"])
        res[[i]] <- res[[i]][keep, ]
    }
    return(res)
}


main <- function() {
    target1 <- paste0(Sys.getenv("rdata_dir"), "txi_", cohort, ".Rda")
    target2 <- paste0(Sys.getenv("ann_dir"), "pheno_", cohort, ".tsv")
    target3 <- paste0(Sys.getenv("genmod_dir"), "herv_ids.tsv")
    inputs <- load_inputs(target1, target2, target3)
    dds <- create_dds(inputs)
    if (cohort != "B") {
        svaobj <- prep_sva(dds)
        dds <- run_sva(svaobj, dds)
    }
    adj <- count_adjustment(dds)
    ddsobj <- setNames(list(dds, adj), c("dds", "adj"))
    target <- paste0(Sys.getenv("rdata_dir"), "dds_", cohort, ".Rda")
    check_dir(Sys.getenv("rdata_dir"))
    write_dds(ddsobj, target)
    ddsobj <- filter_genes(ddsobj$dds, ddsobj$adj, inputs$hervs)
    dgeobj <- test_dge(ddsobj$dds, ddsobj$adj)
    target <- paste0(Sys.getenv("rdata_dir"), "dge_", cohort, ".Rda")
    check_dir(Sys.getenv("rdata_dir"))
    write_dge(dgeobj, target)
    check_dir(Sys.getenv("table_dir"))
    write_gene_list(dgeobj$res)
    res <- extract_hervs(dgeobj$res)
    write_gene_list(res, herv=TRUE)
    print(lapply(dgeobj$res, summary))
}


main()
