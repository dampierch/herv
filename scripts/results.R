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


args <- commandArgs(trailingOnly=TRUE)
cohort <- args[2]
cat(paste("Cohort:", cohort, "\n"))


ggp_theme_default <- theme(
    panel.background=element_rect(fill="white"),
    panel.grid.major=element_line(color="white"),
    panel.grid.minor=element_line(color="white"),
    plot.background=element_rect(fill="white"),
    plot.margin=margin(t=1, r=1, b=1, l=1, unit="lines"),
    plot.title=element_text(size=12, face="bold", hjust=0.5),
    plot.subtitle=element_text(size=10, face="plain", hjust=0.5),
    axis.title.x=element_text(size=10, face="plain"),
    axis.title.y=element_text(size=10, face="plain"),
    axis.text.x=element_text(size=9, face="plain"),
    axis.text.y=element_text(size=9, face="plain"),
    # axis.ticks.x=element_none(),
    # axis.ticks.y=element_none(),
    axis.line.x.bottom=element_line(),
    axis.line.y.left=element_line(),
    legend.key=element_rect(fill="white"),
    legend.position="top",
    legend.title=element_blank(),
    strip.background=element_rect(fill="black"),
    strip.text=element_text(colour="white")
)


ggp_theme_col <- ggp_theme_default +
    theme(
        plot.margin=margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="lines"),
        plot.title=element_text(size=12, face="plain", hjust=0.5),
        axis.text.x=element_text(size=10, face="plain", angle=45, hjust=1),
        axis.text.y=element_text(size=10, face="plain"),
        axis.line.x.bottom=element_blank(),
        axis.line.y.left=element_blank(),
        legend.position="none"
    )


ggp_theme_vol <- ggp_theme_default +
    theme(
        plot.margin=margin(t=0.25, r=0.25, b=0.25, l=0.25, unit="lines"),
        legend.key=element_blank(),
        legend.key.size=unit(0.001, "cm"),
        legend.key.width=unit(0.001, "cm"),
        legend.position=c(0.18, 0.9),
        legend.spacing.x=unit(0.001, "cm"),
        legend.spacing.y=unit(0.001, "cm"),
        legend.text=element_text(size=8, margin=margin(t=0.001))
    )


ggp_theme_box <- ggp_theme_default +
    theme(
        plot.margin=margin(t=0.5, r=1, b=0.5, l=1, unit="lines"),
        axis.text.x=element_text(size=9, face="plain", angle=45, hjust=1),
        axis.text.y=element_text(size=10, face="plain"),
        legend.position="none",
        legend.key=element_blank(),
        legend.key.size=unit(0.75, "cm"),
        legend.key.width=unit(0.75, "cm"),
        legend.spacing.x=unit(0.5, "cm"),
        legend.spacing.y=unit(0.5, "cm"),
        legend.text=element_text(size=9, margin=margin(t=0.1))
    )


draw_key_default <- function(data, params, size) {
    ## this function is here in case GeomBar$draw_key must be manipulated and
    ## then reset once occasion for altered parameters done
    ## https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2
    if (is.null(data$size)) {
        data$size <- 0.5
    }
    lwd <- min(data$size, min(size) / 4)
    grid::rectGrob(
        width = grid::unit(1, "npc") - grid::unit(lwd, "mm"),
        height = grid::unit(1, "npc") - grid::unit(lwd, "mm"),
        gp = grid::gpar(
            col = data$colour %||% NA,
            fill = alpha(data$fill %||% "grey20", data$alpha),
            lty = data$linetype %||% 1,
            lwd = lwd * .pt,
            linejoin = params$linejoin %||% "mitre",
            lineend = if (identical(params$linejoin, "round")) "round" else "square"
        )
    )
}


draw_key_large <- function(data, params, size) {
    ## this function is here in case GeomBar$draw_key must be reset to allow for
    ## larger spacing between legend keys; not used in this script
    ## https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2
    lwd <- min(data$size, min(size) / 4)
    grid::rectGrob(
        width = grid::unit(0.6, "npc"),
        height = grid::unit(0.6, "npc"),
        gp = grid::gpar(
            col = data$colour,
            fill = alpha(data$fill, data$alpha),
            lty = data$linetype,
            lwd = lwd * .pt,
            linejoin = "mitre"
        )
    )
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


load_res_tables <- function(study, scope="all") {
    ## scope can be "all" or "hrv"
    ## order genes by gene_id instead of rank so that order is consistent across
    ## different results tables to facilitate comparisons across tables
    cat("Loading results tables\n")
    if (study == "herv") {
        in_path <- Sys.getenv("table_dir")
        if (cohort == "B") {
            in_name <- c(
                paste0("res_B_CRCvNAT_", scope, ".tsv")
            )
        } else {
            in_name <- c(
                paste0("res_", cohort, "_NATvHLT_", scope, ".tsv"),
                paste0("res_", cohort, "_CRCvHLT_", scope, ".tsv"),
                paste0("res_", cohort, "_CRCvNAT_", scope, ".tsv")
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


plot_cor <- function(x, y, label, r) {
    df <- data.frame(herv=x, field=y)
    ggp_title <- "Pearson Correlation"
    ggp_subtitle <- label
    ggp_xlab <- "L2FC from HERV Study"
    ggp_ylab <- "L2FC from Field Effect Study"
    txt <- paste("R =", round(r, 2))
    ggp <- ggplot(df, aes(x=herv, y=field)) +
        geom_point(size=0.6, stroke=0, shape=16, alpha=0.5) +
        labs(
            title=element_blank(),
            subtitle=ggp_subtitle,
            x=ggp_xlab,
            y=ggp_ylab
        ) +
        annotate("text", x=(min(df$herv) + 1), y=0.95*max(df$field), label=txt,
            size=2.5) +
        ggp_theme_default
    return(ggp)
}


calc_cor <- function(l1, l2, method="pearson") {
    cat("Checking correlation with prior study\n")
    r <- list()
    p <- list()
    for (i in seq_len(length(l1))) {
        intsx <- intersect(l1[[i]]$gene_id, l2[[i]]$gene_id)
        x <- l1[[i]][l1[[i]]$gene_id %in% intsx, "log2FoldChange"]
        y <- l2[[i]][l2[[i]]$gene_id %in% intsx, "log2FoldChange"]
        label <- sub("CRC", "TUM", unlist(strsplit(names(l1)[i], "_"))[3])
        r[[label]] <- cor(x, y, method=method)
        p[[label]] <- plot_cor(x, y, label, r[[label]])
    }
    return(setNames(list(r, p), c("cor", "fig")))
}


write_cor <- function(plotlist, target) {
    p_width <- 9
    p_height <- 3
    fig <- cowplot::plot_grid(
        plotlist[["NATvHLT"]], plotlist[["TUMvNAT"]], plotlist[["TUMvHLT"]],
        nrow=1, ncol=3,
        labels=c("A", "B", "C")
    )
    pdf(file=target, width=p_width, height=p_height)
    print(fig)
    dev.off()
    cat("figure written to", target, "\n")
}


correlation_analysis <- function() {
    tab_herv <- load_res_tables("herv")
    tab_field <- load_res_tables("field")
    l <- calc_cor(tab_herv, tab_field)
    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-field-cor_A.pdf")
    write_cor(l$fig, target)
}


count_herv_ids <- function(gtf, ids) {
    cat("Counting total herv transcripts\n")
    cmd = paste(
        paste("zcat", gtf),
        "awk -v FS=\"\\t\" '!($0~/^#/) && $3 == \"transcript\" {print $0}'",
        "wc -l",
        sep=" | "
    )
    x1 <- as.numeric(system(cmd, intern=TRUE))
    cat("Counting select herv ids\n")
    df <- data.frame(readr::read_tsv(ids))
    x2 <- nrow(df)
    x3 <- length(unique(df$herv_id))
    nms <- c("Total Transcripts", "Select Transcripts", "HERV Genes")
    return(setNames(list(x1, x2, x3), nms))
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


count_herv_cats <- function(cv) {
    nms <- base::unlist(
        base::lapply(
            base::lapply(
                strsplit(cv, "[|]"), unpack_id
            ),
            head, 1
        ),
        use.names=FALSE
    )
    cats <- c(
        "HERVE", "HERVH", "HERVK", "HERVL"
    )
    cnts <- list()
    for (e in cats) {
        cntr <- 0
        for (i in nms) {
            if (i == e) {
                cntr <- cntr + 1
            }
        }
        cnts[[e]] <- cntr
    }
    cnts[["Other"]] <- length(nms) - sum(unlist(cnts))
    return(cnts)
}


fill_params <- function(l, params) {
    if ("Total Transcripts" %in% names(l)) {
        params[["ggp_title"]] <- "Number of HERV Elements Analyzed"
        n <- format(l[["Total Transcripts"]], big.mark=",")
        params[["ggp_sub"]] <- paste("From", n, "starting transcripts")
        l <- l[names(l) != "Total Transcripts"]
    } else {
        params[["ggp_title"]] <- "Classes of HERV Genes Tested"
        params[["ggp_sub"]] <- element_blank()
    }
    params[["ggp_ylab"]] <- "Counts"
    return(setNames(list(l, params), c("l", "params")))
}


format_data <- function(l) {
    df <- data.frame(base::t(data.frame(l, check.names=FALSE)))
    df[, 2] <- rownames(df)
    colnames(df) <- c("counts", "categories")
    if ("Other" %in% rownames(df)) {
        tmp <- df[df$categories != "Other", ]
        df$categories <- factor(
            df$categories,
            levels=c(
                rownames(tmp[order(tmp$counts, decreasing=TRUE), ]),
                "Other"
            )
        )
    } else {
        df$categories <- factor(
            df$categories,
            levels=rownames(df[order(df$counts, decreasing=TRUE), ])
        )
    }
    return(df)
}


plot_hrvsum <- function(inlist) {
    ggp_l <- list()
    cntr <- 1
    for (l in inlist) {
        params <- list()
        int_l <- fill_params(l, params)
        df <- format_data(int_l$l)
        p <- int_l$params
        ggp_l[[cntr]] <- ggplot(df, aes(x=categories, y=counts)) +
            geom_col() +
            labs(title=p[["ggp_title"]], subtitle=p[["ggp_sub"]],
                x=element_blank(), y=p[["ggp_ylab"]]
            ) +
            geom_text(
                stat='identity',
                aes(label=format(counts, big.mark=",")),
                vjust=-0.5,
                size=3
            ) +
            scale_y_continuous(limits=c(0, max(df$counts) * 1.05)) +
            ggp_theme_col
        cntr <- cntr + 1
    }
    return(ggp_l)
}


write_hrvsum <- function(plotlist, target) {
    p_width <- 7
    p_height <- 3.5
    pdf(file=target, width=p_width, height=p_height)
    print(
        cowplot::plot_grid(
            plotlist=plotlist,
            nrow=1,
            ncol=2,
            labels=c("A", "B")
        )
    )
    dev.off()
    cat("figure written to", target, "\n")
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


shrink_effect <- function(dgeobj) {
    ## apparently cannot use apeglm due to third contrast not being present in
    ## resultsNames(dds) and apeglm not taking contrast argument; therefore,
    ## prefer ashr with res object provided instead of coef or contrast argument
    ##
    ## incidentally, the shrinkage has no apparent effect on the plot
    cat("Shrinking effect sizes for visualization\n")
    cons <- list()
    m <- combn(levels(colData(dgeobj$dds)$phenotype), 2)
    for (i in seq_len(ncol(m))) {
        cons[[i]] <- c("phenotype", rev(m[, i]))
        names(cons) <- c(
            names(cons)[seq_len(length(cons) - 1)],
            paste(rev(m[, i]), collapse="v")
        )
    }
    dfs <- list()
    for (i in seq_len(length(cons))) {
        res <- dgeobj$res[[i]][rownames(dgeobj$dds), ]
        tmp <- lfcShrink(
            dgeobj$dds,
            res=res,
            type="ashr"
        )
        if (all(res$gene_id == tmp$gene_id)) {
            res$group <- "Other"
            idx <- res$padj < 0.05 & grepl("HERV", res$gene_id) & tmp$log2FoldChange < 0
            res[idx, "group"] <- "HERV Down"
            idx <- res$padj < 0.05 & grepl("HERV", res$gene_id) & tmp$log2FoldChange > 0
            res[idx, "group"] <- "HERV Up"
            dfs[[i]] <- data.frame(
                gene_id=res$gene_id,
                gene_name=res$gene_name,
                herv=grepl("HERV", res$gene_id),
                effect=tmp$log2FoldChange,
                pval=res$pvalue,
                logp=-log10(res$pvalue),
                padj=res$padj,
                group=res$group
            )
        } else {
            cat("Gene IDs are not equal\n")
        }
    }
    names(dfs) <- names(cons)
    return(dfs)
}


plot_volcano <- function(df, label) {
    ## repeated calls to geom_point ensure HERV point is last point plotted
    ggp_subtitle <- label
    ggp_xlab <- expression("Log"["2"]*" Fold Change")
    ggp_ylab <- expression("-log"["10"]*"("*italic("P")*")")
    ggp <- ggplot(
            data=df,
            aes(x=effect, y=logp, colour=group, size=group, alpha=group)
        ) +
        geom_point(stroke=0, shape=16) +
        # geom_point(size=1, stroke=0, shape=16, alpha=0.5) +
        # geom_point(data=base::subset(df, herv==TRUE & effect>0),
        #     size=1, stroke=0, shape=16, colour="red", alpha=0.25) +
        # geom_point(data=base::subset(df, herv==TRUE & effect<0),
        #     size=1, stroke=0, shape=16, colour="blue", alpha=0.25) +
        geom_point(data=base::subset(df, herv==TRUE & effect>0 & padj<0.05),
            size=1.75, stroke=0, shape=16, colour="red", alpha=0.25) +
        geom_point(data=base::subset(df, herv==TRUE & effect<0 & padj<0.05),
            size=1.75, stroke=0, shape=16, colour="blue", alpha=0.25) +
        labs(
            title=element_blank(),
            subtitle=ggp_subtitle,
            x=ggp_xlab,
            y=ggp_ylab
        ) +
        ggp_theme_vol
    if (cohort %in% c("A", "C")) {
        ggp <- ggp +
            scale_colour_manual(values=c("blue", "red", "grey")) +
            scale_size_manual(values=c(1.75, 1.75, 1)) +
            scale_alpha_manual(values=c(0.5, 0.5, 0.75)) +
            guides(colour=guide_legend(label.position="right",
                override.aes=list(alpha=1, size=1.75))
            )
        if (label == "NATvHLT") {
            ggp <- ggp +
                scale_x_continuous(limits=c(-7, 7)) +
                scale_y_continuous(limits=c(-1, 60))
        } else if (label == "TUMvHLT") {
            ggp <- ggp +
                scale_x_continuous(limits=c(-11, 11)) +
                scale_y_continuous(limits=c(-1, 300))
        } else {
            ggp <- ggp +
                scale_x_continuous(limits=c(-10, 10)) +
                scale_y_continuous(limits=c(-1, 175))
        }
    } else if (cohort == "B") {
        ggp <- ggp +
            scale_colour_manual(values=c("red", "grey")) +
            scale_size_manual(values=c(1.75, 1)) +
            scale_alpha_manual(values=c(0.5, 0.75)) +
            scale_x_continuous(limits=c(-10, 10)) +
            scale_y_continuous(limits=c(-1, 60)) +
            guides(colour=guide_legend(label.position="right",
                override.aes=list(alpha=1, size=1.75))
            )
    }
    return(ggp)
}


assemble_volcano <-function(dgeobj) {
    cat("Assembling effect size vs significance (i.e. volcano) plots\n")
    dfs <- shrink_effect(dgeobj)
    p <- list()
    for (i in seq_len(length(dfs))) {
        label <- sub("CRC", "TUM", names(dfs)[i])
        p[[label]] <- plot_volcano(dfs[[i]], label)
    }
    return(p)
}


write_volcano <- function(plotlist, target) {
    if (length(plotlist) > 1) {
        p_width <- 9
        p_height <- 3
        fig <- cowplot::plot_grid(
            plotlist[["NATvHLT"]], plotlist[["TUMvNAT"]], plotlist[["TUMvHLT"]],
            nrow=1,
            ncol=3,
            labels=c("A", "B", "C")
        )
    } else {
        p_width <- 3
        p_height <- 3
        fig <- cowplot::plot_grid(
            plotlist=plotlist,
            nrow=1, ncol=1, labels=c("A")
        )
    }
    pdf(file=target, width=p_width, height=p_height)
    print(fig)
    dev.off()
    cat("figure written to", target, "\n")
}


volcano_analysis <- function(dgeobj) {
    pl <- assemble_volcano(dgeobj)
    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-volcano_", cohort, ".pdf")
    write_volcano(pl, target)
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


map_cytoband <- function(tophervs) {
    target <- paste0(Sys.getenv("ref_dir"), "annotation/ucsc/",
        "cytoBand.txt.gz"
    )
    con <- gzfile(target, "rb")
    cytomap <- readr::read_tsv(con,
        col_names=c("chrom", "start", "stop", "band", "giemsa")
    )
    close(con)
    l <- setNames(base::lapply(strsplit(tophervs, "[|]"), unpack_id), tophervs)
    tophervs2 <- base::lapply(l, parse_band, cytomap)
    return(tophervs2)
}


select_tophervs <- function(model="naive") {
    ## model indicates underlying assumption about gene expression in NAT vs HLT
    ## field or naive; field assumes NAT is intermediate between HLT and TUM;
    ## naive assumes NAT == HLT
    cat("Selecting top hits for box plots\n")
    restabs <- load_res_tables("herv", scope="hrv")
    nms <- character()
    for (i in seq_len(length(restabs))) {
        nms <- c(nms,
            sub("CRC", "TUM", unlist(strsplit(names(restabs)[i], "_"))[3])
        )
    }
    names(restabs) <- nms
    conds <- list()
    conds[[1]] <- restabs[["NATvHLT"]]$padj < 0.05
    conds[[2]] <- restabs[["NATvHLT"]]$log2FoldChange > 0
    conds[[3]] <- restabs[["TUMvNAT"]]$padj < 0.05
    conds[[4]] <- restabs[["TUMvNAT"]]$log2FoldChange > 0
    conds[[5]] <- restabs[["TUMvHLT"]]$padj < 0.05
    conds[[6]] <- restabs[["TUMvHLT"]]$log2FoldChange > 0
    if (model == "field") {
        idx <- conds[[1]] & conds[[2]] & conds[[3]] & conds[[4]]  ## A: 1 gene
    } else {
        idx <- conds[[3]] & conds[[4]] & conds[[5]] & conds[[6]]  ## A: 11 genes
    }
    tophervs <- restabs[["TUMvHLT"]][idx, "gene_id"]
    tophervs2 <- map_cytoband(tophervs)
    return(tophervs2)
}


format_counts <- function(tophervs, dgeobj) {
    cat("Formatting counts for box plots\n")
    idx <- rownames(dgeobj$adj) %in% names(tophervs)
    mx <- dgeobj$adj[idx, ]
    dimnames(mx)[[1]] <- unlist(tophervs, use.names=FALSE)
    dimnames(mx)[[2]] <- gsub("CRC", "TUM",
        as.character(colData(dgeobj$dds)$phenotype)
    )
    df <- reshape2::melt(mx)
    colnames(df) <- c("herv_id", "phenotype", "count")
    df$phenotype <- factor(df$phenotype, levels=c("HLT", "NAT", "TUM"))
    dat <- dplyr::group_by(base::subset(df, phenotype=="TUM"), herv_id) %>%
        dplyr::summarise(mean=mean(count)) %>% dplyr::arrange(dplyr::desc(mean))
    ord <- as.character(dat$herv_id)
    df$herv_id <- factor(df$herv_id, levels=ord)
    if (cohort == "A") {
        df1 <- df[df$herv_id %in% ord[1:6], ]
        df1$herv_id <- factor(df1$herv_id)
        df2 <- df[df$herv_id %in% ord[7:11], ]
        df2$herv_id <- factor(df2$herv_id)
        dfs <- setNames(list(df1, df2), c("upper", "lower"))
    } else {
        dfs <- list(df)
    }
    return(dfs)
}


plot_counts <- function(dfs) {
    cat("Making box plots\n")
    color_values <- c("blue", "grey", "red")
    ggp_title <- "Expression Levels of Select HERV Genes"
    ggp_ylab <- "Scaled Counts"
    lidx <- 0
    pl <- list()
    for (df in dfs) {
        lidx <- lidx + 1
        ggp <- ggplot(df, aes(x=herv_id, y=count, colour=phenotype)) +
            geom_boxplot(outlier.size=-1, position=position_dodge2()) +
            geom_point(size=0.25, alpha=0.5, position=position_jitterdodge()) +
            labs(title=ggp_title, x=element_blank(), y=ggp_ylab) +
            scale_colour_manual(values=color_values) +
            ggp_theme_box
        if (cohort == "A" & length(levels(df$herv_id)) < 6) {
            ggp <- ggp +
                theme(legend.position="right") +
                labs(title=element_blank())
        }
        pl[[lidx]] <- ggp
    }
    return(pl)
}


write_counts <- function(plotlist, target) {
    if (length(plotlist) > 1) {
        p_width <- 7
        p_height <- 5
        fig <- cowplot::plot_grid(
            plotlist=plotlist,
            nrow=2,
            ncol=1,
            labels=c("A", NULL)
        )
    } else {
        p_width <- 5
        p_height <- 3
        fig <- cowplot::plot_grid(
            plotlist=plotlist,
            nrow=1,
            ncol=1,
            labels=c("A")
        )
    }
    pdf(file=target, width=p_width, height=p_height)
    print(fig)
    dev.off()
    cat("figure written to", target, "\n")
}


tophervs_analysis <- function(dgeobj) {
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
        tophervs_analysis(dgeobj)
    }
}


main()




# target <- "/scratch/chd5n/test-fig.pdf"
# ggsave(target, p, device="pdf", width=15, height=4, units="in")
