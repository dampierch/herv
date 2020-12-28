## to be sourced in results.R


write_cor_select <- function(plotlist, target) {
    p_width <- 6
    p_height <- 3
    fig <- cowplot::plot_grid(
        plotlist[["TUMvHLT"]], plotlist[["TUMvNAT"]],
        nrow=1, ncol=2,
        labels=c("A", "B")
    )
    pdf(file=target, width=p_width, height=p_height)
    print(fig)
    dev.off()
    cat("figure written to", target, "\n")
    return(fig)
}


fig_supp_cor_field <- function() {
    cat("Making supp fig HERV vs field effect corr\n")
    tab_herv <- load_res_tables("herv")
    tab_field <- load_res_tables("field")
    l <- calc_cor(tab_herv, tab_field)
    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-field-cor-select_A.pdf")
    fig <- write_cor_select(l$fig, target)
    return(fig)
}


load_inputs_prelim_expr <- function(target1, target2, target3) {
    cat("Loading inputs for prelim expr\n")
    df <- data.frame(ensembldb::genes(EnsDb.Hsapiens.v86))
    idx <- df$gene_biotype == "protein_coding"
    pcg <- as.vector(df[idx, "gene_id"])
    hv <- data.frame(readr::read_tsv(target1))
    v <- as.vector(hv$herv_id)
    load(target2, verbose=TRUE)
    load(target3, verbose=TRUE)
    return(setNames(
        list(pcg, v, ddsobj, dgeobj),
        c("pcg", "hervs", "ddsobj", "dgeobj")
    ))
}


plot_prelim_expr <- function(pcg, hervs, ddsobj, pheno=NULL) {
    cat("Plotting prelim expr\n")
    ggp_title <- "Protein Coding and HERV Genes"
    ggp_xlab <- "Expression Level (All)"
    ggp_ylab <- "Number of Genes"
    idx <- rownames(ddsobj$adj) %in% pcg
    adjpcg <- ddsobj$adj[idx, ]
    idx <- rownames(ddsobj$adj) %in% hervs
    adjhrv <- ddsobj$adj[idx, ]
    if (!is.null(pheno)) {
        ggp_xlab <- "Expression Level (TUM)"
        idx <- colData(ddsobj$dds)$phenotype == pheno
        adjpcg <- adjpcg[ , idx]
        adjhrv <- adjhrv[ , idx]
    }
    df <- rbind(
        data.frame(expr=apply(adjpcg, 1, median), lab="GENCODE"),
        data.frame(expr=apply(adjhrv, 1, median), lab="HERV")
    )
    ggp <- ggplot(df, aes(x=expr, fill=lab)) +
        geom_histogram(
            position="identity",
            binwidth=0.6,
            alpha=0.5
        ) +
        geom_vline(
            xintercept=median(subset(df, lab="GENCODE")$expr),
            linetype="dashed", size=0.25
        ) +
        geom_vline(
            xintercept=diff(c(
                sd(subset(df, lab="GENCODE")$expr),
                median(subset(df, lab="GENCODE")$expr)
            )),
            linetype="dotted", size=0.25
        ) +
        labs(
            subtitle=ggp_title,
            x=ggp_xlab,
            y=ggp_ylab
        ) +
        scale_fill_grey(
            start=0.8,
            end=0.2
        ) +
        ggp_theme_default +
        theme(plot.margin=margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="lines"))
    return(ggp)
}


run_hrvnums2 <- function(hervs, ddsobj, dgeobj) {
    cat("Running the numbers for HERVs\n")
    x0 <- 1001931
    x1 <- length(hervs)
    x2 <- length(unique(hervs))
    x3 <- sum(rownames(ddsobj$adj) %in% hervs)
    x4 <- sum(rownames(dgeobj$adj) %in% hervs)
    labs <- c("Total Transcripts", "Selected Transcripts",
        "Selected HERV Genes", "Detected HERV Genes", "Tested HERV Genes")
    return(setNames(list(x0, x1, x2, x3, x4), labs))
}


write_prelim_expr <- function(plotlist, target) {
    p_width <- 6
    p_height <- 6
    fig <- cowplot::plot_grid(
        plotlist=plotlist,
        nrow=2,
        ncol=2,
        labels=c("A", "B", "C", "D")
    )
    pdf(file=target, width=p_width, height=p_height)
    print(fig)
    dev.off()
    cat("figure written to", target, "\n")
    return(fig)
}


fig_prelim_expr <- function() {
    cat("Making fig of HERVs studied\n")
    target1 <- paste0(Sys.getenv("genmod_dir"), "herv_ids.tsv")
    target2 <- paste0(Sys.getenv("rdata_dir"), "dds_A.Rda")
    target3 <- paste0(Sys.getenv("rdata_dir"), "dge_A.Rda")
    l <- load_inputs_prelim_expr(target1, target2, target3)
    hrvnums <- run_hrvnums2(l$hervs, l$ddsobj, l$dgeobj)
    hrvnames <- rownames(l$dgeobj$adj)[grepl("HERV", rownames(l$dgeobj$adj))]
    hrvcats <- count_herv_cats(hrvnames)
    sumlist <- setNames(list(hrvnums, hrvcats), c("HERV Genes", "HERV Classes"))
    pl <- plot_hrvsum(sumlist)
    pl[[1]] <- pl[[1]] +
        labs(title="HERV Elements Studied") +
        coord_cartesian(clip='off')
    pl[[3]] <- plot_prelim_expr(l$pcg, l$hervs, l$ddsobj)
    pl[[4]] <- plot_prelim_expr(l$pcg, l$hervs, l$ddsobj, pheno="CRC")
    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-prelim-expr_A.pdf")
    fig <- write_prelim_expr(pl, target)
    return(fig)
}


plot_volcano_select <- function(plotlist, target) {
    cwp <- cowplot::plot_grid(
        plotlist[["TUMvHLT"]], plotlist[["TUMvNAT"]],
        nrow=1,
        ncol=2,
        labels=c("A", "B")
    )
    return(cwp)
}


plot_counts_select <- function(plotlist, target) {
    cwp <- cowplot::plot_grid(
        plotlist=plotlist,
        nrow=2,
        ncol=1,
        labels=c("C", NULL)
    )
    return(cwp)
}


write_dge_result <- function(plotlist, target) {
    p_width <- 6
    p_height <- 9
    fig <- cowplot::plot_grid(
        plotlist=plotlist,
        nrow=2,
        ncol=1,
        labels=c(NULL, NULL),
        rel_heights=c(1, 2)
    )
    pdf(file=target, width=p_width, height=p_height)
    print(fig)
    dev.off()
    cat("figure written to", target, "\n")
}


fig_dge_result <- function() {
    cat("Making fig of HERV DGE results\n")
    cwpl <- list()

    target <- paste0(Sys.getenv("rdata_dir"), "dge_A.Rda")
    load(target, verbose=TRUE)
    pl <- assemble_volcano(dgeobj)
    pl[["TUMvHLT"]] <- pl[["TUMvHLT"]] + labs(subtitle="TUM - HLT")
    pl[["TUMvNAT"]] <- pl[["TUMvNAT"]] + labs(subtitle="TUM - NAT")
    cwpl[[1]] <- plot_volcano_select(pl)

    tophervs <- select_tophervs()
    dfs <- format_counts(tophervs, dgeobj)
    pl <- plot_counts(dfs)
    pl[[1]] <- pl[[1]] +
        labs(title=element_blank(),
            subtitle="Expression Levels of Tumor-Specific HERV Genes"
        )
    cwpl[[2]] <- plot_counts_select(pl)

    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-dge-result.pdf")
    fig <- write_dge_result(cwpl, target)
    return(fig)
}
