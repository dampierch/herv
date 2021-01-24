## to be sourced in results.R


read_charac_inputs <- function(cohort) {
    target <- paste0(Sys.getenv("rdata_dir"), "dge_", cohort, ".Rda")
    dgeobj <- load_dgeobj(target)
    restabs <- fill_restabs()
    l <- extract_main_hervs(restabs)
    herv_ids <- intersect(
        l[["A"]][["res_A_CRCvHLT_hrv"]], l[["A"]][["res_A_CRCvNAT_hrv"]]
    )
    names(herv_ids) <- map_cytoband(herv_ids)
    names(herv_ids) <- gsub("HERV", "", names(herv_ids))
    names(herv_ids) <- gsub("_", " ", names(herv_ids))
    ## CTLA4 not among genes tested
    # markers <- c("ENSG00000120217", "ENSG00000188389")
    # names(markers) <- c("CD274", "PDCD1")
    markers <- c(
        "ENSG00000145649", "ENSG00000100453", "ENSG00000100450",
        "ENSG00000113088", "ENSG00000197540", "ENSG00000180644"
    )
    names(markers) <- c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1")
    return(setNames(
        list(herv_ids, dgeobj, markers),
        c("herv_ids", "dgeobj", "markers")
    ))
}


prep_charac_df <- function(herv_ids, dgeobj, markers, feature) {
    i <- rownames(dgeobj$adj) %in% c(herv_ids)
    j <- colData(dgeobj$dds)$phenotype == "CRC"
    if (feature %in% c("msi_status", "stage")) {
        j1 <- !is.na(colData(dgeobj$dds)[ , feature])
        j <- j & j1
    } else if (feature == "immunomarker") {
        i <- rownames(dgeobj$adj) %in% c(herv_ids, markers)
    } else if (feature == "age") {
        j <- !is.na(colData(dgeobj$dds)$age_days)
    } else {
        warning("feature unrecognized, default data.frame provided")
    }
    df <- data.frame(dgeobj$adj[i, j])
    colnames(df) <- colData(dgeobj$dds)$subject_id[j]
    return(df)
}


explore_charac_markers <- function(herv_ids, markers, df) {
    ## draft code to determine what to show
    ## for command line manual execution, not to run
    l_cor <- list()
    for (hid in herv_ids) {
        for (mrk in markers) {
            id <- paste(
                names(herv_ids)[hid == herv_ids],
                names(markers)[mrk == markers],
                sep="::"
            )
            l_cor[[id]] <- cor(unlist(df[hid, ]), unlist(df[mrk, ]))
        }
    }
    a <- unlist(l_cor)
    df1 <- data.frame(id=names(a), cor=as.numeric(a))
    df1$herv <- unlist(apply(df1, 1, fun <- function(r) {
        unlist(strsplit(r[1], "::"))[1]
    }))
    df1$marker <- unlist(apply(df1, 1, fun <- function(r) {
        unlist(strsplit(r[1], "::"))[2]
    }))
    df1[abs(df1$cor) > 0.25, ]
    df1[order(df1$cor), ]
}


prep_charac_im1 <- function(herv_ids, markers, df) {
    cd8 <- apply(df[markers, ], 2, mean)
    df <- df[herv_ids, ]
    df$herv_id <- names(herv_ids[match(rownames(df), herv_ids)])
    df <- reshape2::melt(
        df, id.vars="herv_id", variable.name="subject_id", value.name="herv"
    )
    df$marker <- cd8[match(df$subject_id, names(cd8))]
    df$group <- df$herv_id
    return(df)
}


prep_charac_im2 <- function(herv_ids, markers, df, hid) {
    hrv <- unlist(df[herv_ids[hid], ])
    df <- df[markers, ]
    df$gene_name <- names(markers[match(rownames(df), markers)])
    df <- reshape2::melt(
        df, id.vars="gene_name", variable.name="subject_id", value.name="marker"
    )
    df$herv <- hrv[match(df$subject_id, names(hrv))]
    df$group <- df$gene_name
    return(df)
}


plot_charac_im <- function(df, herv_levels, hid=NULL) {
    if ("gene_name" %in% colnames(df)) {
        ggp_title <- paste("HERV", hid, "and Cytotoxic Activity")
        nrow <- 2
        ggp_theme_charac <- ggp_theme_charac +
            theme(strip.text=element_text(colour="white", face="italic"))
    } else {
        df$group <- factor(df$group, levels=herv_levels)
        ggp_title <- "HERVs and Cytotoxic Activity Index"
        nrow <- 4
    }
    ggp_ylab <- "Marker Expression"
    ggp_xlab <- "HERV Expression"
    ggp <- ggplot(df, aes(x=herv, y=marker, group=group)) +
        geom_point(size=0.25, alpha=0.5) +
        geom_smooth(method="lm", se=FALSE, color=muted("red"), size=0.8) +
        labs(subtitle=ggp_title, x=ggp_xlab, y=ggp_ylab) +
        facet_wrap(vars(group), nrow=nrow, scales="free") +
        ggp_theme_charac
    return(ggp)
}


make_charac_im_plotlist <- function(l, df, herv_levels) {
    pl <- list()
    df1 <- prep_charac_im1(l$herv_ids, l$markers, df)
    hid <- "1 I 3p22.3"
    df2 <- prep_charac_im2(l$herv_ids, l$markers, df, hid)
    pl[[1]] <- plot_charac_im(df1, herv_levels)
    pl[[2]] <- plot_charac_im(df2, herv_levels, hid=hid)
    return(pl)
}


write_charac_im <- function(target, pl) {
    pdf(target, height=8, width=4)
    print(
        cowplot::plot_grid(
            plotlist=pl, labels=c("a", "b"), nrow=2, rel_heights=c(1, 0.5)
        )
    )
    dev.off()
    cat("plot written to", target, "\n")
}


prep_charac_msi_stage_age <- function(herv_ids, dgeobj, df, feature) {
    df$herv_id <- names(herv_ids[match(rownames(df), herv_ids)])
    df <- reshape2::melt(
        df, id.vars="herv_id", variable.name="subject_id", value.name="count"
    )
    idx <- match(df$subject_id, colData(dgeobj$dds)$subject_id)
    if (feature == "msi_status") {
        df$msi_status <- colData(dgeobj$dds)[idx, "msi_status"]
        df <- df[df$msi_status != "Indeterminate", ]
        df$x <- df$msi_status
    } else if (feature == "stage") {
        df$stage <- toupper(colData(dgeobj$dds)[idx, "stage"])
        df$x <- df$stage
    }  else {
        df$age <- colData(dgeobj$dds)[idx, "age_days"]/365.25
        df$phenotype <- colData(dgeobj$dds)[idx, "phenotype"]
    }
    return(df)
}


plot_charac_msi_stage <- function(df, herv_levels) {
    df$herv_id <- factor(df$herv_id, levels=herv_levels)
    if ("msi_status" %in% colnames(df)) {
        ggp_title <- "HERVs and MSI Status"
    } else {
        ggp_title <- "HERVs and Tumor Stage"
    }
    ggp_ylab <- "Expression Level"
    ggp <- ggplot(df, aes(x=x, y=count, colour=x)) +
        geom_boxplot(outlier.size=-1, position=position_dodge2()) +
        geom_point(size=0.25, alpha=0.5, position=position_jitterdodge()) +
        labs(subtitle=ggp_title, x=element_blank(), y=ggp_ylab) +
        facet_wrap(vars(herv_id), nrow=4, scales="free_y") +
        scale_colour_grey(start=0.8, end=0.2) +
        ggp_theme_charac +
        theme(axis.text.x=element_text(size=9, face="plain", angle=45, hjust=1))
    return(ggp)
}


write_charac_msi_stage <- function(target, pl) {
    pdf(target, height=6, width=8)
    print(
        cowplot::plot_grid(
            plotlist=pl, labels=c("a", "b"), nrow=1
        )
    )
    dev.off()
    cat("plot written to", target, "\n")
}


plot_charac_age_corr <- function(df, herv_levels) {
    df$herv_id <- factor(df$herv_id, levels=herv_levels)
    ggp_title <- "HERVs and Age"
    ggp_ylab <- "Expression Level"
    ggp_xlab <- "Years"
    ggp <- ggplot(df, aes(x=age, y=count)) +
        geom_point(size=0.25, alpha=0.5) +
        geom_smooth(method="lm", se=FALSE, color=muted("red"), size=0.8) +
        labs(subtitle=ggp_title, x=ggp_xlab, y=ggp_ylab) +
        facet_wrap(vars(herv_id), nrow=4, scales="free") +
        ggp_theme_charac
    target <- "/scratch/chd5n/test_age.pdf"
    ggsave(target, ggp, height=6, width=4, unit="in")
    return(ggp)
}


prep_charac_age_lm <- function(df) {
    cor_labs <- list()
    fits <- list()
    resids <- list()
    for (each in unique(df$herv_id)) {
        d <- df[df$herv_id == each, ]
        cor_labs[[each]] <- cor(d$count, d$age)
        df_coeff <- data.frame(
            summary(lm(count ~ age + phenotype, data=d))$coefficients
        )
        colnames(df_coeff) <- c("beta", "se", "tval", "pval")
        rownames(df_coeff) <- c("int", "AGE", "NAT", "TUM")
        df_coeff$predictor <- c("int", "AGE", "NAT", "TUM")
        fits[[each]] <- df_coeff
        df_resid <- data.frame(
            resid=lm(count ~ age + phenotype, data=d)$residuals,
            herv_id=each
        )
        resids[[each]] <- df_resid
    }
    return(setNames(list(fits, resids), c("fits", "resids")))
}


check_charac_residuals_freqpoly <- function(df) {
    ## this function takes a dataframe of residuals from a lm fit and generates
    ## frequency distributions for selected genes
    ggp_title <- "Distribution of Residuals"
    ggp_xlab <- "Residuals"
    ggp_ylab <- "Frequency"
    bwd <- 0.8
    ggp <- ggplot(df,aes(x=resid, colour=herv_id)) +
        geom_freqpoly(binwidth=bwd, na.rm=TRUE) +
        labs(title=ggp_title, x=ggp_xlab, y=ggp_ylab) +
        ggp_theme_default +
        theme(legend.position="none")
    target <- "/scratch/chd5n/test_pred.pdf"
    ggsave(target, ggp, device="pdf", height=3, width=3, unit="in")
    cat("plot saved to", target ,"\n")
}


plot_charac_age_tval <- function(df) {
    ## this function takes a dataframe of summary stats from a lm fit and
    ## generates box plots of test stats for selected genes across predictors
    ggp_title <- "HERVs and LM Predictors"
    ggp_xlab <- "Predictors"
    ggp_ylab <- "Test Statistic"
    ggp <- ggplot(df, aes(x=predictor, y=tval)) +
        geom_boxplot(outlier.size=-1, width=0.5) +
        geom_jitter(size=1, width=0.2) +
        geom_hline(yintercept=0, colour=muted("red"), linetype="dashed") +
        labs(subtitle=ggp_title, x=ggp_xlab, y=ggp_ylab) +
        ggp_theme_box
    cwp <- cowplot::plot_grid(ggp, NULL, nrow=2, labels=c(NULL, NULL))
    return(cwp)
}


write_charac_age <- function(target, pl) {
    pdf(target, height=6, width=6)
    print(
        cowplot::plot_grid(
            plotlist=pl, labels=c("a", "b"), nrow=1, rel_widths=c(0.66, 0.33)
        )
    )
    dev.off()
    cat("plot written to", target, "\n")
}
