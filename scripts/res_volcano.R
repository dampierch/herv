## to be sourced in results.R


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
