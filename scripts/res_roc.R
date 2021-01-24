library(readr)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(pROC)


source("res_utilities.R")
source("res_themes.R")
source("res_val.R")


args <- commandArgs(trailingOnly=TRUE)
cohort <- args[2]
cat(paste("Cohort:", cohort, "\n"))


read_inputs <- function(cohort) {
    target <- paste0(Sys.getenv("rdata_dir"), "dge_", cohort, ".Rda")
    dgeobj <- load_dgeobj(target)
    target <- paste0(Sys.getenv("rdata_dir"), "txi_", cohort, ".Rda")
    load(target, verbose=TRUE)
    restabs <- fill_restabs()
    return(setNames(list(dgeobj, txi, restabs), c("dgeobj", "txi", "restabs")))
}


set_stage <- function(truth) {
    truth[is.na(truth$stage_r), "stage_r"] <- "nan"
    truth[truth$phenotype %in% c("HLT", "NAT"), "stage_r"] <- "hlt"
    l <- list()
    l[["nan"]] <- NA
    l[["hlt"]] <- 0
    l[["i"]] <- 1
    l[["ii"]] <- 2
    l[["iii"]] <- 3
    l[["iii/iv"]] <- 4
    l[["iv"]] <- 4
    truth$stage <- unlist(
        lapply(truth$stage_r, fun <- function(x) {l[[x]]}),
        use.names=FALSE
    )
    return(truth)
}


set_truth <- function(dgeobj) {
    truth <- data.frame(
        subject_id=colData(dgeobj$dds)$subject_id,
        phenotype=colData(dgeobj$dds)$phenotype,
        sex=colData(dgeobj$dds)$sex,
        age=colData(dgeobj$dds)$age_days/365.25,
        stage_r=colData(dgeobj$dds)$stage
    )
    truth <- set_stage(truth)
    return(truth)
}


set_predictor <-function(dgeobj, txi, restabs) {
    l <- extract_main_hervs(restabs)
    herv_ids <- intersect(
        l[["A"]][["res_A_CRCvHLT_hrv"]], l[["A"]][["res_A_CRCvNAT_hrv"]]
    )
    tpms <- txi$abundance[match(herv_ids, rownames(txi$abundance)), ]
    colnames(tpms) <- colData(dgeobj$dds)$subject_id
    return(tpms)
}


fit_classifier <- function(df, tpms, g1, g2="CRC") {
    ## need these herv_id variables due to eval(parse(text=hid)) code in glm
    herv_ids1 <- rownames(tpms)
    ## usual id
    names(herv_ids1) <- map_cytoband(herv_ids1)
    names(herv_ids1) <- gsub("HERV", "", names(herv_ids1))
    names(herv_ids1) <- gsub("_", " ", names(herv_ids1))
    ## glm acceptable id
    rownames(tpms) <- map_cytoband(rownames(tpms))
    rownames(tpms) <- gsub(" ", ".", rownames(tpms))
    herv_ids2 <- rownames(tpms)
    names(herv_ids2) <- names(herv_ids1)
    ## set truth values for binary classifier
    truth <- df[as.character(df$phenotype) %in% c(g1, g2), ]
    truth$phenotype <- as.numeric(factor(truth$phenotype)) - 1
    ## set dataframe for glm
    data <- cbind(truth, t(tpms[ , as.character(df$phenotype) %in% c(g1, g2)]))
    fit <- list()
    # roc <- list()  ## if want to rank by auc
    for (hid in herv_ids2) {
        fit[[hid]] <- glm(
            phenotype ~ eval(parse(text=hid)),
            data=data, family=binomial
        )
        ## roc will give auc for ordering by auc if desired
        # roc[[hid]] <- pROC::roc(
        #     data$phenotype, fit[[hid]]$fitted.values, plot=FALSE
        # )
    }
    ## switch names back to usual style
    names(fit) <- names(herv_ids2[match(herv_ids2, names(fit))])
    # o <- order(
    #     unlist(lapply(roc, fun <- function(x) {as.numeric(x[["auc"]])})),
    #     decreasing=TRUE
    # )
    # fit1 <- fit[o]
    return(setNames(list(truth, fit), c("truth", "fit")))
}


print_roc <- function(target, truth, fit) {
    ## this function is for draft exploration of data
    ## uses base r graphic parameters
    pdf(target, height=8, width=6)
    par(mfrow=c(4, 3), mai=c(0.25, 0.25, 0.5, 0.25))
    par(pty="s")
    lapply(fit, fun <- function(x) {
        pROC::roc(
            truth$phenotype, x$fitted.values, plot=TRUE, legacy.axes=TRUE,
            xlab="False Positive Rate", ylab="True Positive Rate",
            main=x$herv_id, font.main=1,
            col="#377eb8", lwd=2,
            print.auc=TRUE, print.auc.cex=0.9
        )
    })
    dev.off()
    cat("plot written to", target, "\n")
}


plot_roc <- function(truth, fit) {
    ## make roc objects
    roc_list <- lapply(fit, fun <- function(x) {
        pROC::roc(
            truth$phenotype, x$fitted.values, plot=FALSE, legacy.axes=TRUE
        )
    })
    ## extract coordinates
    dfs <- lapply(roc_list, fun <- function(x) {
        df <- pROC::coords(x, "all", transpose=FALSE)
        return(df[rev(seq(nrow(df))), ])
    })
    ## extract aucs
    aucs <- lapply(roc_list, fun <- function(x) {as.numeric(x$auc)})
  	## add herv_ids
  	for (n in names(dfs)) {dfs[[n]]$herv_id <- n}
  	## combine all
  	df <- do.call(rbind, dfs)
    ## set levels
    herv_levels <- c(
        "H Xp22.32", "IP10F 9q13", "H 13q33.3", "K9 13q12.11", "IP10F 9p11.2",
        "H 20p11.23", "1 I 3p22.3", "L 9q21.11", "H 5q15", "E a 7q22.1",
        "H 13q14.11"
    )
    df$herv_id <- factor(df$herv_id, levels=herv_levels)
    ## set up auc labels
    a <- dplyr::bind_rows(aucs)
    df_lab <- data.frame(herv_id=names(a), auc=as.numeric(a))
    df_lab$herv_id <- factor(df_lab$herv_id, levels=herv_levels)
    ## make plot object
    ggp <- ggplot(df, aes(x=1-specificity, y=sensitivity, group=herv_id)) +
        geom_line(size=1, colour="#377eb8", alpha=0.75) +
        geom_abline(intercept=0, slope=1, size=0.2, linetype="dashed") +
        geom_text(
            data=df_lab,
            mapping=aes(
                x=0.75, y=0.4,
                label=paste("AUC", format(round(auc, 3), nsmall=3))
            ),
            size=3
        ) +
        labs(x="False Positive Rate", y="True Positive Rate") +
        scale_x_continuous(
            limits=c(0, 1), labels=fun <- function(x) {round(x, 1)}
        ) +
        scale_y_continuous(
            limits=c(0, 1), labels=fun <- function(x) {round(x, 1)}
        ) +
        facet_wrap(vars(herv_id), nrow=4) + ggp_theme_charac
    return(ggp)
}


write_roc <- function(target, pl, labs=NULL, wd=4.5) {
    if (length(pl) > 1) {
        labs <- c("a", "b")
        wd <- 9
    }
    pdf(target, height=6, width=wd)
    print(
        cowplot::plot_grid(
            plotlist=pl, labels=labs, nrow=1
        )
    )
    dev.off()
    cat("plot written to", target, "\n")
}


main <- function(cohort) {
    cat("Running ROC analysis for cohort", cohort, "\n")
    inputs <- read_inputs(cohort)
    truth <- set_truth(inputs$dgeobj)
    tpms <- set_predictor(inputs$dgeobj, inputs$txi, inputs$restabs)
    if (cohort %in% c("A", "C")) {
        groups <- c("HLT", "NAT")
    } else {
        groups <- c("NAT")
    }
    pl <- list()
    for (g1 in groups) {
        l <- fit_classifier(truth, tpms, g1)
        pl[[g1]] <- plot_roc(l$truth, l$fit)
    }
    target <- paste0(Sys.getenv("plot_dir"), "herv-roc_", cohort, ".pdf")
    write_roc(target, pl)
}


main(cohort)
