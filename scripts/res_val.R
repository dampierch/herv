## to be sourced in results.R


fill_restabs <- function() {
    cat("Loading all results for validation testing\n")
    restabs <- list()
    for (i in c("A", "B", "C")) {
        restabs[[i]] <- load_res_tables("herv", scope="hrv", set_cohort=i)
    }
    return(restabs)
}


calc_cor_val <- function(restabs, cohort1, cohort2, method="pearson") {
    cat("Checking correlations between HERV cohorts\n")
    s1 <- paste0("res_", cohort1, "_CRCvNAT_hrv")
    s2 <- paste0("res_", cohort2, "_CRCvNAT_hrv")
    intx <- intersect(
        restabs[[cohort1]][[s1]]$gene_id,
        restabs[[cohort2]][[s2]]$gene_id
    )
    i <- restabs[[cohort1]][[s1]]$gene_id %in% intx
    lfc1 <- restabs[[cohort1]][[s1]][i, "log2FoldChange"]
    i <- restabs[[cohort2]][[s2]]$gene_id %in% intx
    lfc2 <- restabs[[cohort2]][[s2]][i, "log2FoldChange"]
    cc <- cor(lfc1, lfc2, method=method)
    return(setNames(list(lfc1, lfc2, cc), c("lfc1", "lfc2", "coef")))
}


plot_cor_val <- function(l, labs) {
    df <- data.frame(x=l$lfc1, y=l$lfc2)
    ggp_title <- "Pearson Correlation"
    ggp_subtitle <- paste("HERVs Tested =", length(l$lfc1))
    ggp_xlab <- paste("L2FC from", labs[1])
    ggp_ylab <- paste("L2FC from", labs[2])
    txt <- paste("R =", round(l$coef, 2))
    ggp <- ggplot(df, aes(x, y)) +
        geom_point(size=1.5, stroke=0, shape=16, alpha=0.5) +
        labs(
            title=element_blank(),
            subtitle=ggp_subtitle,
            x=ggp_xlab,
            y=ggp_ylab
        ) +
        annotate("text", x=(min(df$x) + 1), y=0.95*max(df$y), label=txt,
            size=2.5) +
        ggp_theme_default
    return(ggp)
}


val_corr_analysis <- function(restabs) {
    cat("Performing correlations tests for validation\n")
    l <- calc_cor_val(restabs, "A", "B")
    labs <- c("Discovery", "Ind Cohort #2")
    ggpAB <- plot_cor_val(l, labs)

    l <- calc_cor_val(restabs, "C", "B")
    labs <- c("Ind Cohort #1", "Ind Cohort #2")
    ggpCB <- plot_cor_val(l, labs)

    l <- calc_cor_val(restabs, "A", "C")
    labs <- c("Discovery", "Ind Cohort #1")
    ggpAC <- plot_cor_val(l, labs)

    pl <- list(ggpAC, ggpAB, ggpCB)
    return(pl)
}


write_cor_val <- function(plotlist, target) {
    cat("Writing HERV cohort correlation plots\n")
    p_width <- 9
    p_height <- 3
    fig <- cowplot::plot_grid(
        plotlist=plotlist,
        nrow=1, ncol=3,
        labels=c("A", "B", "C")
    )
    pdf(file=target, width=p_width, height=p_height)
    print(fig)
    dev.off()
    cat("figure written to", target, "\n")
}


extract_main_hervs <- function(restabs) {
    cat("Extracting HERV names for main results\n")
    idx <- lapply(restabs, fun <- function(x1) {
        idx1 <- lapply(x1, fun <- function(x2) {
            a <- x2$padj < 0.05
            b <- x2$log2FoldChange > 0
            c <- a & b
            return(c)
        })
        return(idx1)
    })
    l <- list()
    for (n1 in names(restabs)) {
        l[[n1]] <- list()
        for (n2 in names(restabs[[n1]])) {
            l[[n1]][[n2]] <- restabs[[n1]][[n2]][idx[[n1]][[n2]], "gene_id"]
        }
    }
    return(l)
}


stringent_val <- function(l) {
    cat("Performing stringent validation testing\n")
    topA <- intersect(
        l[["A"]][["res_A_CRCvHLT_hrv"]], l[["A"]][["res_A_CRCvNAT_hrv"]]
    )
    topB <- l[["B"]][["res_B_CRCvNAT_hrv"]]
    topC <- intersect(
        l[["C"]][["res_C_CRCvHLT_hrv"]], l[["C"]][["res_C_CRCvNAT_hrv"]]
    )
    val1 <- intersect(topA, topB)
    val2 <- intersect(topA, topC)
    return(setNames(list(val1, val2), c("ind2", "ind1")))
}


build_val_df <- function(restabs, l) {
    cat("Building dataframe of validated results\n")
    topA <- intersect(
        l[["A"]][["res_A_CRCvHLT_hrv"]], l[["A"]][["res_A_CRCvNAT_hrv"]]
    )
    valtabs <- lapply(restabs, fun <- function(x1) {
        vt <- lapply(x1, fun <- function(x2) {
            a <- x2$gene_id %in% topA
            b <- x2$pvalue < 0.05
            c <- x2[a & b, ]
            return(c)
        })
        return(vt)
    })
    valdf <- rbind(valtabs$C$res_C_CRCvHLT_hrv, valtabs$C$res_C_CRCvNAT_hrv)
    fields <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
    valdf <- valdf[!duplicated(valdf$gene_id), fields]
    valdf$fc <- round(2 ^ valdf$log2FoldChange)
    valdf$CTL <- unlist(
        lapply(valdf$gene_id, fun <- function(x) {
            a <- c()
            if (x %in% valtabs$C$res_C_CRCvHLT_hrv$gene_id) {a <- append(a, "HLT")}
            if (x %in% valtabs$C$res_C_CRCvNAT_hrv$gene_id) {a <- append(a, "NAT")}
            a <- paste(a, collapse=" & ")
            return(a)
        }), use.names=FALSE
    )
    valdf$gene_cyt <- unlist(map_cytoband(valdf$gene_id), use.names=FALSE)
    valdf$order <- order(
        as.numeric(
            unlist(
                lapply(
                    strsplit(
                        gsub("q", ".6",
                            gsub("p", ".5",
                                gsub("\\.", "_",
                                    gsub("X", 23,
                                        unlist(
                                            lapply(strsplit(valdf$gene_cyt, " "), "[", 2),
                                            use.names=FALSE
                                        )
                                    )
                                )
                            )
                        ), "_"
                    ), "[", 1
                ), use.names=FALSE
            )
        )
    )
    valdf <- valdf[valdf$order, ]
    fields <- c("HERV ID", "Mean", "Log2(FC)", "SE", "Stat", "P", "P Adj", "FC", "CTL", "HERV", "o")
    colnames(valdf) <- fields
    fields <- c("HERV", "Mean", "FC", "Log2(FC)", "SE", "Stat", "P", "P Adj", "CTL")
    valdf <- valdf[ , fields]
    for (i in c("Mean", "Log2(FC)", "SE", "Stat", "P", "P Adj")) {
        valdf[ , i] <- signif(valdf[ , i], 2)
    }
    valdf$HERV <- gsub("HERV", "", valdf$HERV)
    valdf$HERV <- gsub("_", " ", valdf$HERV)
    return(valdf)
}


write_valtable <- function(valdf, target) {
    cat("Writing validated results to tex table\n")
    obj <- xtable::xtable(
        valdf,
        caption="Tumor-specific HERV genes validated in independent cohort",
        label="tab_val",
        digits=c(0, 0, 0, 0, 2, 2, 2, 2, 2, 0),
        display=c("d", "s", "f", "f", "f", "f", "f", "E", "E", "s")
    )
    xtable::print.xtable(obj, file=target, size="small", include.rownames=FALSE,
        hline.after=c(0,nrow(obj)), caption.placement="top",
        sanitize.colnames.function=function(x){paste('{\\textbf{',x,'}}', sep='')}
    )
    cat("table written to", target, "\n")
}
