## to be sourced in results.R


plot_cor_select <- function(plotlist) {
    cwp <- cowplot::plot_grid(
        plotlist[["TUMvHLT"]], plotlist[["TUMvNAT"]],
        nrow=1, ncol=2,
        labels=c("A", "B")
    )
    return(cwp)
}


plot_cor_val_select <- function(plotlist) {
    cwp <- cowplot::plot_grid(
        plotlist[[1]], plotlist[[2]],
        nrow=1, ncol=2,
        labels=c("C", "D")
    )
    return(cwp)
}


write_cor_select <- function(plotlist, target) {
    p_width <- 6
    p_height <- 6
    fig <- cowplot::plot_grid(
        plotlist=plotlist,
        nrow=2, ncol=1,
        labels=c(NULL, NULL)
    )
    pdf(file=target, width=p_width, height=p_height)
    print(fig)
    dev.off()
    cat("figure written to", target, "\n")
    return(fig)
}


fig_supp_cor_field_val <- function() {
    cat("Making supp fig HERV corrs\n")
    cwpl <- list()

    tab_herv <- load_res_tables("herv")
    tab_field <- load_res_tables("field")
    l <- calc_cor(tab_herv, tab_field)
    l$fig[["TUMvHLT"]] <- l$fig[["TUMvHLT"]] + labs(subtitle="TUM - HLT")
    l$fig[["TUMvNAT"]] <- l$fig[["TUMvNAT"]] + labs(subtitle="TUM - NAT")
    cwpl[[1]] <- plot_cor_select(l$fig)

    restabs <- fill_restabs()
    pl <- val_corr_analysis(restabs)
    cwpl[[2]] <- plot_cor_val_select(pl)

    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-cor-select.pdf")
    fig <- write_cor_select(cwpl, target)
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
    ggp_title <- "GENCODE Genes and HERV Loci"
    ggp_xlab <- "Expression Level (All Samples)"
    ggp_ylab <- "Number of Genes/Loci"
    idx <- rownames(ddsobj$adj) %in% pcg
    adjpcg <- ddsobj$adj[idx, ]
    idx <- rownames(ddsobj$adj) %in% hervs
    adjhrv <- ddsobj$adj[idx, ]
    if (!is.null(pheno)) {
        ggp_xlab <- "Expression Level (Tumor Samples)"
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
        "Selected HERV Loci", "Detected HERV Loci", "Tested HERV Loci")
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
    sumlist <- setNames(list(hrvnums, hrvcats), c("HERV Loci", "HERV Groups"))
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
    return(fig)
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
            subtitle="Expression Levels of Tumor Specific HERV Loci"
        )
    cwpl[[2]] <- plot_counts_select(pl)

    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-dge-result.pdf")
    fig <- write_dge_result(cwpl, target)
    return(fig)
}


get_herv_cat <- function(hervids) {
    maincats <- c("HERVE", "HERVH", "HERVK", "HERVL")
    maincls <- c("E", "H", "K", "L")
    nms <- unname(sapply(lapply(strsplit(hervids, "[|]"), unpack_id), head, 1))
    a <- sapply(nms, fun <- function(x) {
        if (x %in% maincats) {
            y <- plyr::mapvalues(x, from=maincats, to=maincls, warn_missing=FALSE)
        } else {
            y <- "Other"
        }
        return(y)
    }, USE.NAMES=FALSE)
    return(a)
}


get_heat_data <- function(target) {
    cat("Getting data ready for heatmap\n")
    load(target, verbose=TRUE)

    hrvadj <- dgeobj$adj[grepl("HERV", rownames(dgeobj$adj)), ]
    hrvmn <- apply(hrvadj, 1, mean)
    hrvsd <- apply(hrvadj, 1, sd)
    hrvz <- data.frame((hrvadj - hrvmn) / hrvsd)
    colnames(hrvz) <- colData(dgeobj$dds)$subject_id
    hrvz <- hrvz[ , order(colData(dgeobj$dds)$phenotype)]
    hrvz$herv_id <- rownames(hrvz)
    hrvz$herv_cl <- get_herv_cat(hrvz$herv_id)
    hrvz <- hrvz[order(hrvz$herv_cl, decreasing=TRUE), ]
    hrvz$n <- 0
    for (i in levels(as.factor(hrvz$herv_cl))) {
        hrvz[hrvz$herv_cl == i, "n"] <- seq_len(nrow(hrvz[hrvz$herv_cl == i, ]))
    }
    hrvz$y <- paste(hrvz$herv_cl, hrvz$n, sep="_")
    hrvz <- hrvz[ , !(colnames(hrvz) == "n")]

    df <- melt(
        hrvz,
        id.vars=c("herv_id", "herv_cl", "y"),
        variable.name="subj_id",
        value.name="Z"
    )
    df$pheno <- plyr::mapvalues(
          df$subj_id,
          from=colData(dgeobj$dds)$subject_id,
          to=as.character(colData(dgeobj$dds)$phenotype)
    )
    df$pheno <- as.character(df$pheno)
    df[df$pheno == "CRC", "pheno"] <- "TUM"
    df$n <- 0
    for (i in levels(as.factor(df$pheno))) {
        df[df$pheno == i, "n"] <- plyr::mapvalues(
            as.character(df[df$pheno == i, "subj_id"]),
            from=unique(as.character(df[df$pheno == i, "subj_id"])),
            to=seq_len(length(unique(as.character(df[df$pheno == i, "subj_id"]))))
        )
    }
    df$x <- paste(df$pheno, df$n, sep="_")
    df <- df[ , !(colnames(df) == "n")]
    return(df)
}


plot_heat_big <- function(df) {
    cat("Making heatmap\n")
    ggp <- ggplot(df, aes(x=x, y=y, fill=Z)) +
        geom_raster() +
        scale_fill_gradient(low="white", high="black", name="Z-score") +
        scale_y_discrete(limits=rev(levels(as.factor(df$y)))) +
        coord_cartesian(clip="off") +
        ggp_theme_heat +
        theme(
            axis.text.x=element_blank(),
            axis.text.y=element_blank()
        )

    N <- list()
    text <- list()
    par <- list()
    N$HLT <- max(as.numeric(unlist(lapply(strsplit(df[df$pheno=="HLT", "x"], "_"), "[", 2))))
    N$NAT <- max(as.numeric(unlist(lapply(strsplit(df[df$pheno=="NAT", "x"], "_"), "[", 2))))
    N$TUM <- max(as.numeric(unlist(lapply(strsplit(df[df$pheno=="TUM", "x"], "_"), "[", 2))))
    text$HLT <- grid::textGrob("HLT", gp=grid::gpar(fontsize=10))
    text$NAT <- grid::textGrob("NAT", gp=grid::gpar(fontsize=10))
    text$TUM <- grid::textGrob("TUM", gp=grid::gpar(fontsize=10))
    par$HLT <- c((N$HLT / 2) - 5, (N$HLT / 2) + 5, -5, -5)
    par$NAT <- c(N$HLT + (N$NAT / 2) - 5, N$HLT + (N$NAT / 2) + 5, -5, -5)
    par$TUM <- c(N$HLT + N$NAT + (N$TUM / 2) - 5, N$HLT + N$NAT + (N$TUM / 2) + 5,
        -5, -5)
    ggp <- ggp +
        annotation_custom(text$HLT, xmin=par$HLT[1], xmax=par$HLT[2], ymin=par$HLT[3], ymax=par$HLT[4]) +
        annotation_custom(text$NAT, xmin=par$NAT[1], xmax=par$NAT[2], ymin=par$NAT[3], ymax=par$NAT[4]) +
        annotation_custom(text$TUM, xmin=par$TUM[1], xmax=par$TUM[2], ymin=par$TUM[3], ymax=par$TUM[4])

    N1 <- list()
    text1 <- list()
    par1 <- list()
    N1$E <- max(as.numeric(unlist(lapply(strsplit(df[df$herv_cl=="E", "y"], "_"), "[", 2))))
    N1$H <- max(as.numeric(unlist(lapply(strsplit(df[df$herv_cl=="H", "y"], "_"), "[", 2))))
    N1$K <- max(as.numeric(unlist(lapply(strsplit(df[df$herv_cl=="K", "y"], "_"), "[", 2))))
    N1$L <- max(as.numeric(unlist(lapply(strsplit(df[df$herv_cl=="L", "y"], "_"), "[", 2))))
    N1$Other <- max(as.numeric(unlist(lapply(strsplit(df[df$herv_cl=="Other", "y"], "_"), "[", 2))))
    text1$E <- grid::textGrob("E", rot=0, gp=grid::gpar(fontsize=10))
    text1$H <- grid::textGrob("H", rot=0, gp=grid::gpar(fontsize=10))
    text1$K <- grid::textGrob("K", rot=0, gp=grid::gpar(fontsize=10))
    text1$L <- grid::textGrob("L", rot=0, gp=grid::gpar(fontsize=10))
    text1$Other <- grid::textGrob("O", rot=0, gp=grid::gpar(fontsize=10))
    par1$Other <- c(-15, -15, (N1$Other / 2) - 5, (N1$Other / 2) + 5)
    par1$L <- c(-15, -15, N1$Other + (N1$L / 2) - 5, N1$Other + (N1$L / 2) + 5)
    par1$K <- c(-15, -15, N1$Other + N1$L + (N1$K / 2) - 5,
        N1$Other + N1$L + (N1$K / 2) + 5)
    par1$H <- c(-15, -15, N1$Other + N1$L + N1$K + (N1$H / 2) - 5,
        N1$Other + N1$L + N1$K + (N1$H / 2) + 5)
    par1$E <- c(-15, -15, N1$Other + N1$L + N1$K + N1$H + (N1$E / 2) - 5,
        N1$Other + N1$L + N1$K + N1$H + (N1$E / 2) + 5)
    ggp <- ggp +
        annotation_custom(text1$E, xmin=par1$E[1], xmax=par1$E[2], ymin=par1$E[3], ymax=par1$E[4]) +
        annotation_custom(text1$H, xmin=par1$H[1], xmax=par1$H[2], ymin=par1$H[3], ymax=par1$H[4]) +
        annotation_custom(text1$K, xmin=par1$K[1], xmax=par1$K[2], ymin=par1$K[3], ymax=par1$K[4]) +
        annotation_custom(text1$L, xmin=par1$L[1], xmax=par1$L[2], ymin=par1$L[3], ymax=par1$L[4]) +
        annotation_custom(text1$Other, xmin=par1$Other[1], xmax=par1$Other[2],
            ymin=par1$Other[3], ymax=par1$Other[4])
    return(ggp)
}


plot_heat_small <- function(df) {
    ## this turned out to be too problematic to interpret and is deprectated
    df <- df %>% dplyr::group_by(herv_cl, pheno) %>%
        dplyr::summarise(Z=median(Z), n=n())
    ggp <- ggplot(df, aes(x=pheno, y=herv_cl, fill=Z)) +
        geom_tile() +
        scale_fill_gradient(low="white", high="black", name="Z-score") +
        scale_y_discrete(limits=rev(levels(as.factor(df$herv_cl)))) +
        ggp_theme_heat
    return(ggp)
}


write_heat <- function(plotlist, target) {
    p_width <- 6   ## 7 for A,B
    p_height <- 3  ## 3.5 for A,B
    fig <- cowplot::plot_grid(
        plotlist=plotlist,
        nrow=1,
        ncol=1,
        # labels=c("A", "B"),
        labels=c(NULL)
        # rel_widths=c(1, 0.6)
    )
    pdf(file=target, width=p_width, height=p_height)
    print(fig)
    dev.off()
    cat("figure written to", target, "\n")
    return(fig)
}


fig_heat <- function() {
    pl <- list()
    target <- paste0(Sys.getenv("rdata_dir"), "dge_A.Rda")
    df <- get_heat_data(target)
    pl[[1]] <- plot_heat_big(df)
    pl[[2]] <- plot_heat_small(df)
    target_dir <- Sys.getenv("plot_dir")
    check_dir(target_dir)
    target <- paste0(target_dir, "herv-heat.pdf")
    fig <- write_heat(pl, target)
    return(fig)
}


devs_ordering <- function() {
    ## miscellaneous code that may be helpful some day for ordering data
    idx <- list()
    for (i in levels(as.factor(df$pheno))) {
        idx[[i]] <- df$pheno == i
    }

    start_num <- 0
    for (i in seq_len(length(idx))) {
        df[idx[[i]], "xorder"] <- seq.int(
            max(c(start_num, max(df$xorder))) + 1,
            max(c(start_num, max(df$xorder))) + sum(idx[[i]])
        )
    }

    idx <- list()
    for (i in levels(as.factor(df$herv_cl))) {
        idx[[i]] <- df$herv_cl == i
    }

    start_num <- 0
    for (i in seq_len(length(idx))) {
        df[idx[[i]], "yorder"] <- seq.int(
            max(c(start_num, max(df$yorder))) + 1,
            max(c(start_num, max(df$yorder))) + sum(idx[[i]])
        )
    }
}
