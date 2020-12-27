## to be sourced in results.R


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


set_paths <- function(in_path, in_name) {
    targets <- list()
    for (e in in_name) {
        target <- paste0(in_path, e)
        n <- substr(e, 1, (nchar(e) - 4))
        if (file.exists(target)) {
            targets[[n]] <- target
        } else {
            cat("problem locating file", target, "\n")
        }
    }
    return(targets)
}


load_res_tables <- function(study, scope="all", set_cohort=NULL) {
    ## scope can be "all" or "hrv"
    ## order genes by gene_id instead of rank so that order is consistent across
    ## different results tables to facilitate comparisons across tables
    cat("Loading results tables\n")
    if (is.null(set_cohort)) {
        cohort_choice <- cohort
    } else {
        cohort_choice <- set_cohort
    }
    if (study == "herv") {
        in_path <- Sys.getenv("table_dir")
        if (cohort_choice == "B") {
            in_name <- c(
                paste0("res_B_CRCvNAT_", scope, ".tsv")
            )
        } else {
            in_name <- c(
                paste0("res_", cohort_choice, "_NATvHLT_", scope, ".tsv"),
                paste0("res_", cohort_choice, "_CRCvHLT_", scope, ".tsv"),
                paste0("res_", cohort_choice, "_CRCvNAT_", scope, ".tsv")
            )
        }
    } else {
        in_path <- "/nv/t74/cphgdesk/share/cphg_CaseyLab/Papers/Field_Effect/tables/"
        in_name <- c(
            "p-nat-hlt-ResTab.csv",
            "p-crc-hlt-ResTab.csv",
            "p-crc-nat-ResTab.csv"
        )
    }
    targets <- set_paths(in_path, in_name)
    tables <- list()
    for (i in names(targets)) {
        target <- targets[[i]]
        if (study == "herv") {
            df <- data.frame(readr::read_tsv(target))
            df <- df[order(df$gene_id), ]
        } else {
            df <- data.frame(readr::read_csv(target))
            colnames(df) <- c("gene_id", colnames(df)[2:ncol(df)])
            df <- df[order(df$gene_id), ]
        }
        tables[[i]] <- df
    }
    return(tables)
}


load_dgeobj <- function(target) {
    ## dgeobj is a list with dds, adj, res objects
    ## dds is deseq2 data set, which includes pheno in colData()
    ## adj is a matrix of vst (and for cohorts A and C batch adjusted) counts
    ## res is a deseq2 results object
    cat("Loading dgeobj\n")
    nm <- base::load(target, verbose=TRUE)
    return(dgeobj)
}


unpack_id <- function(e) {
    ## unpacks information in herv_id
    x <- unlist(strsplit(e[2], "-"))[1]
    hvnm_s <- substr(x, 1, 5)
    hvnm_l <- x
    chrom <- paste0("chr", e[3])
    start <- as.numeric(e[4])
    stop <- as.numeric(e[5])
    return(setNames(c(hvnm_s, hvnm_l, chrom, start, stop),
        c("hrvnm_s", "hrvnm_l", "chrom", "start", "stop"))
    )
}


parse_band <- function(e, cytomap) {
    ## used within map_cytoband()
    cond1 <- cytomap$chrom == e["chrom"]
    cond2 <- cytomap$start <= as.numeric(e["start"])
    cond3 <- cytomap$stop >= as.numeric(e["start"])
    idx <- cond1 & cond2 & cond3
    return(paste(
        e["hrvnm_l"],
        paste0(strsplit(e["chrom"], "chr")[[1]][2], cytomap[idx, "band"])
    ))
}


map_cytoband <- function(hervlist) {
    target <- paste0(Sys.getenv("ref_dir"), "annotation/ucsc/",
        "cytoBand.txt.gz"
    )
    con <- gzfile(target, "rb")
    cytomap <- readr::read_tsv(con,
        col_names=c("chrom", "start", "stop", "band", "giemsa")
    )
    close(con)
    l <- setNames(base::lapply(strsplit(hervlist, "[|]"), unpack_id), hervlist)
    hervlist2 <- base::lapply(l, parse_band, cytomap)
    return(hervlist2)
}
