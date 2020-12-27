## to be sourced in results.R


simple_qq <- function(restabs, target) {
    cat("Making simple QQ\n")
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


bacon_qq <- function(tab_herv, target=NULL) {
    cat("Running BACON\n")
    bc <- lapply(tab_herv, fun <- function(x) {bacon(
            effectsizes=x$log2FoldChange,
            standarderrors=x$lfcSE,
    )})
    if (!is.null(target)) {
        cat("Making BACON QQ\n")
        pdf(file=target, height=3, width=6)
        print(lapply(bc, fun <- function(x) {
            plot(x, type="hist")
        }))
        print(lapply(bc, fun <- function(x) {
            plot(x, type="qq")
        }))
        dev.off()
        cat("figure written to", target, "\n")
    }
    return(bc)
}


bacon_correction <- function(tab_herv) {
    tab_herv_bac <- lapply(tab_herv, fun <- function(x) {
        pbac <- bacon::pval(
            bacon::bacon(
                effectsizes=x$log2FoldChange,
                standarderrors=x$lfcSE,
            )
        )
        x$pbac <- pbac
        return(x)
    })
    return(tab_herv_bac)
}


bonferroni_correction <- function(tab_herv) {
    tab_herv_bon <- lapply(tab_herv, fun <- function(x) {
        pbon <- p.adjust(x$pvalue, method="bonferroni")
        x$pbon <- pbon
        return(x)
    })
    return(tab_herv_bon)
}


check_results <- function(tab_herv, type) {
    if (type == "basic") {
        l <- tab_herv
        p <- "padj"
    } else if (type == "bacon") {
        l <- bacon_correction(tab_herv)
        p <- "pbac"
    } else if (type == "bonferroni") {
        l <- bonferroni_correction(tab_herv)
        p <- "pbon"
    }
    c1 <- l$res_A_CRCvHLT_all[ , p] < 0.05
    c2 <- grepl("HERV", l$res_A_CRCvHLT_all$gene_id)
    tum_hlt <- l$res_A_CRCvHLT_all[c1 & c2, ]
    c1 <- l$res_A_CRCvNAT_all[ , p] < 0.05
    c2 <- grepl("HERV", l$res_A_CRCvNAT_all$gene_id)
    tum_nat <- l$res_A_CRCvNAT_all[c1 & c2, ]
    sel <- intersect(tum_hlt$gene_id, tum_nat$gene_id)
    c1 <- l$res_A_CRCvNAT_all$gene_id %in% sel
    c2 <- l$res_A_CRCvNAT_all$stat > 0
    print(l$res_A_CRCvNAT_all[c1 & c2, ])
}


check_results2 <- function(tab_herv) {
    bc <- bacon_qq(tab_herv)
    l <- lapply(bc, pval)
    l2 <- lapply(l, fun <- function(p) {p.adjust(p, method="fdr")})
    a <- lapply(l, fun <- function(x) {sum(x<0.05)})
    b <- lapply(l2, fun <- function(x) {sum(x<0.05)})
    cat("simple bacon\n")
    print(a)
    cat("bacon plus fdr\n")
    print(b)
}


screen_results <- function(tab_herv) {
    check_results(tab_herv, "basic")
    check_results(tab_herv, "bacon")
    check_results(tab_herv, "bonferroni")
    check_results2(tab_herv)
}
