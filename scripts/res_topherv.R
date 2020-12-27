## to be sourced in results.R


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
