## to be sourced in results.R


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
