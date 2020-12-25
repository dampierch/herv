## to be sourced in results.R


simple_qq <- function(restabs, target) {
    y <- list()
    y[[1]] <- -log10(restabs[[1]]$pvalue)
    y[[2]] <- -log10(restabs[[2]]$pvalue)
    y[[3]] <- -log10(restabs[[3]]$pvalue)
    x <- -log10(qunif(ppoints(length(y[[1]]))))

    xlab <- expression("Expected -log"["10"]*"("*italic("P")*")")
    ylab <- expression("Observed -log"["10"]*"("*italic("P")*")")

    pdf(file=target, height=3, width=9)
    par(mfrow=c(1,3))
    qqplot(x, y[[1]], xlab=xlab, ylab=ylab, main="NATvHLT", col="blue")
    qqplot(x, y[[3]], xlab=xlab, ylab=ylab, main="TUMvNAT", col="blue")
    qqplot(x, y[[2]], xlab=xlab, ylab=ylab, main="TUMvHLT", col="blue")
    dev.off()
    cat("figure written to", target, "\n")
}


tab_herv <- load_res_tables("herv")
bc <- lapply(tab_herv, fun <- function(x) {bacon(
        effectsizes=x$log2FoldChange,
        standarderrors=x$lfcSE,
)})

target_dir <- Sys.getenv("plot_dir")
check_dir(target_dir)
target <- paste0(target_dir, "herv-qq2_", cohort, ".pdf")
pdf(file=target, height=3, width=6)
print(lapply(bc, fun <- function(x) {
    plot(x, type="hist")
}))
print(lapply(bc, fun <- function(x) {
    plot(x, type="qq")
}))
dev.off()

tab_herv2 <- lapply(tab_herv, fun <- function(x) {
    padj2 <- bacon::pval(
        bacon::bacon(
            effectsizes=x$log2FoldChange,
            standarderrors=x$lfcSE,
        )
    )
    x$padj2 <- padj2
    return(x)
})

lapply(tab_herv2, fun <- function(x) {
    print(head(x[order(x$rank), ]))
})


tab_herv3 <- lapply(tab_herv, fun <- function(x) {
    padj3 <- p.adjust(x$pvalue, method="bonferroni")
    x$padj3 <- padj3
    return(x)
})


c1 <- tab_herv$res_A_CRCvHLT_all$padj < 0.05
c2 <- grepl("HERV", tab_herv$res_A_CRCvHLT_all$gene_id)
tum_hlt <- tab_herv$res_A_CRCvHLT_all[c1 & c2, ]
c1 <- tab_herv$res_A_CRCvNAT_all$padj < 0.05
c2 <- grepl("HERV", tab_herv$res_A_CRCvNAT_all$gene_id)
tum_nat <- tab_herv$res_A_CRCvNAT_all[c1 & c2, ]
sel <- intersect(tum_hlt$gene_id, tum_nat$gene_id)
c1 <- tab_herv$res_A_CRCvNAT_all$gene_id %in% sel
c2 <- tab_herv$res_A_CRCvNAT_all$stat > 0
tab_herv$res_A_CRCvNAT_all[c1 & c2, ]

c1 <- tab_herv2$res_A_CRCvHLT_all$padj2 < 0.05
c2 <- grepl("HERV", tab_herv2$res_A_CRCvHLT_all$gene_id)
tum_hlt <- tab_herv2$res_A_CRCvHLT_all[c1 & c2, ]
c1 <- tab_herv2$res_A_CRCvNAT_all$padj2 < 0.05
c2 <- grepl("HERV", tab_herv2$res_A_CRCvNAT_all$gene_id)
tum_nat <- tab_herv2$res_A_CRCvNAT_all[c1 & c2, ]
sel <- intersect(tum_hlt$gene_id, tum_nat$gene_id)
c1 <- tab_herv2$res_A_CRCvNAT_all$gene_id %in% sel
c2 <- tab_herv2$res_A_CRCvNAT_all$stat > 0
tab_herv2$res_A_CRCvNAT_all[c1 & c2, ]

c1 <- tab_herv3$res_A_CRCvHLT_all$padj3 < 0.05
c2 <- grepl("HERV", tab_herv3$res_A_CRCvHLT_all$gene_id)
tum_hlt <- tab_herv3$res_A_CRCvHLT_all[c1 & c2, ]
c1 <- tab_herv3$res_A_CRCvNAT_all$padj3 < 0.05
c2 <- grepl("HERV", tab_herv3$res_A_CRCvNAT_all$gene_id)
tum_nat <- tab_herv3$res_A_CRCvNAT_all[c1 & c2, ]
sel <- intersect(tum_hlt$gene_id, tum_nat$gene_id)
c1 <- tab_herv3$res_A_CRCvNAT_all$gene_id %in% sel
c2 <- tab_herv3$res_A_CRCvNAT_all$stat > 0
tab_herv3$res_A_CRCvNAT_all[c1 & c2, ]



l <- lapply(bc, pval)
lapply(l, fun <- function(x) {sum(x<0.05)})
l2 <- lapply(l, fun <- function(p) {p.adjust(p, method="fdr", n=length(p))})
