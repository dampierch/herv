library(readr)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(pROC)
library(survival)
library(survminer)


source("res_utilities.R")
source("res_val.R")


read_inputs <- function() {
    ann <- list()
    tpm <- list()
    target <- paste0(Sys.getenv("ann_dir"), "tcga/ann.tsv")
    ann$tcga <- readr::read_tsv(target)
    target <- paste0(Sys.getenv("ann_dir"), "tcga/cbioportal.tsv")
    ann$cbio <- readr::read_tsv(target)
    target <- paste0(Sys.getenv("ann_dir"), "pheno_A.tsv")
    ann$a <- readr::read_tsv(target)
    target <- paste0(Sys.getenv("ann_dir"), "pheno_C.tsv")
    ann$c <- readr::read_tsv(target)
    target <- paste0(Sys.getenv("ann_dir"), "pheno_B.tsv")
    ann$b <- readr::read_tsv(target)
    target <- paste0(Sys.getenv("rdata_dir"), "txi_A.Rda")
    load(target, verbose=TRUE)
    tpm$a <- txi$abundance
    target <- paste0(Sys.getenv("rdata_dir"), "txi_C.Rda")
    load(target, verbose=TRUE)
    tpm$c <- txi$abundance
    target <- paste0(Sys.getenv("rdata_dir"), "txi_B.Rda")
    load(target, verbose=TRUE)
    tpm$b <- txi$abundance
    restabs <- fill_restabs()
    return(setNames(list(ann, tpm, restabs), c("ann", "tpm", "restabs")))
}


parse_inputs <- function(ann, tpm, restabs) {
    ## extract ids of hervs of interest
    l <- extract_main_hervs(restabs)
    herv_ids <- intersect(
        l[["A"]][["res_A_CRCvHLT_hrv"]], l[["A"]][["res_A_CRCvNAT_hrv"]]
    )
    ## extract tpms of hervs of interest
    tpm <- lapply(tpm, fun <- function(x) {
        x[match(herv_ids, rownames(x)), ]
    })
    ## make herv names easy to read
    tpm <- lapply(tpm, fun <- function(x) {
        rownames(x) <- map_cytoband(herv_ids)
        rownames(x) <- gsub("HERV", "", rownames(x))
        rownames(x) <- gsub("_", " ", rownames(x))
        return(x)
    })
    ## extract ids of tcga samples of interest
    idx <- ann$tcga$tissue_type == "Primary Tumor"
    uuids <- unlist(ann$tcga[idx, "file_id"], use.names=FALSE)
    ## find columns of tcga samples of interest
    idx <- lapply(ann, fun <- function(x) {
        if ("uuid" %in% colnames(x)) {
            match(uuids, x$uuid)[!is.na(match(uuids, x$uuid))]
        }
    })
    ## extract tpms of tcga samples of interest
    tpm1 <- list()
    for (n in names(tpm)) {
        tpm1[[n]] <- tpm[[n]][ , idx[[n]]]
        colnames(tpm1[[n]]) <- unlist(
            ann[[n]][idx[[n]], "subject_id"], use.names=FALSE
        )
    }
    tpms <- cbind(tpm1[[1]], cbind(tpm1[[2]], tpm1[[3]]))
    ## extract vital status of tcga samples of interest
    idx <- match(colnames(tpms), ann$tcga$subject_id)
    fields <- c(
        "subject_id", "vital_status", "days_to_death", "stage", "msi_status",
        "age_at_index", "bmi", "race", "sex"
    )
    truth <- ann$tcga[idx, fields]
    idx <- match(colnames(tpms), unlist(ann$cbio[ , "Patient ID"]))
    fields <- c(
        "Patient ID",
        "Diagnosis Age",
        "Neoplasm Disease Stage American Joint Committee on Cancer Code",
        "Disease Free (Months)",
        "Disease Free Status",
        "Months of disease-specific survival",
        "Disease-specific Survival status",
        "Overall Survival (Months)",
        "Overall Survival Status",
        "Progress Free Survival (Months)",
        "Progression Free Status"
    )
    truth <- cbind(truth, ann$cbio[idx, fields])
    return(setNames(list(tpms, truth), c("tpms", "truth")))
}


fit_classifier <- function(
    truth,
    tpms,
    stage=c("not reported", "stage i", "stage ia", "stage ii", "stage iia",
        "stage iib", "stage iic", "stage iii", "stage iiia", "stage iiib",
        "stage iiic", "stage iv", "stage iva", "stage ivb"
    ))
{
    idx <- truth$stage %in% stage
    truth <- as.numeric(factor(truth$vital_status[idx])) - 1
    pred <- tpms[ , idx]
    fit <- list()
    roc <- list()
    for (hid in rownames(tpms)) {
        fit[[hid]] <- glm(truth ~ pred[hid, ], family=binomial)
        fit[[hid]]$herv_id <- hid
        roc[[hid]] <- pROC::roc(truth, fit[[hid]]$fitted.values, plot=FALSE)
    }
    o <- order(
        unlist(lapply(roc, fun <- function(x) {as.numeric(x[["auc"]])})),
        decreasing=TRUE
    )
    fit1 <- fit[o]
    return(setNames(list(truth, fit1), c("truth", "fit")))
}


print_roc <- function(target, truth, fit) {
    pdf(target, height=8, width=6)
    par(mfrow=c(4, 3), mai=c(0.25, 0.25, 0.5, 0.25))
    par(pty="s")
    lapply(fit, fun <- function(x) {
        pROC::roc(
            truth, x$fitted.values, plot=TRUE, legacy.axes=TRUE,
            xlab="False Positive Rate", ylab="True Positive Rate",
            main=x$herv_id, font.main=1,
            col="#377eb8", lwd=2,
            print.auc=TRUE, print.auc.cex=0.9
        )
    })
    dev.off()
    cat("plot written to", target, "\n")
}


main <- function() {
    l0 <- read_inputs()
    l1 <- parse_inputs(l0$ann, l0$tpm, l0$restabs)
    l2 <- fit_classifier(l1$truth, l1$tpms)
    target <- paste0(
        Sys.getenv("plot_dir"), "herv-surv-roc.pdf"
    )
    print_roc(target, l2$truth, l2$fit)
}


cox_models <- function(l1$truth, l1$tpms) {

    df <- cbind(data.frame(l1$truth), data.frame(t(l1$tpms)))

    preds <- c("sex", "age_at_index", "X1.I.3p22.3", "E.a.7q22.1", "H.5q15",
        "H.13q14.11", "H.20p11.23", "H.13q33.3", "H.Xp22.32", "IP10F.9p11.2",
        "IP10F.9q13", "K9.13q12.11", "L.9q21.11")
    fit <- list()
    for (i in preds) {
        fit[[i]] <- coxph(
            Surv(
                Months.of.disease.specific.survival,
                as.numeric(factor(Disease.specific.Survival.status)) - 1
            ) ~ eval(as.name(i)), data=df
        )
    }

    preds <- c("X1.I.3p22.3", "E.a.7q22.1", "H.5q15",
        "H.13q14.11", "H.20p11.23", "H.13q33.3", "H.Xp22.32", "IP10F.9p11.2",
        "IP10F.9q13", "K9.13q12.11", "L.9q21.11")
    fit <- list()
    for (i in preds) {
        fit[[i]] <- coxph(
            Surv(
                Overall.Survival..Months.,
                as.numeric(factor(Overall.Survival.Status)) - 1
            ) ~ eval(as.name(i)) + age_at_index, data=df
        )
    }

    fivenum(df$H.20p11.23)

    df$group <- "Low"
    df[df$H.20p11.23 > 100, "group"] <- "High"
    fit1 <- coxph(
        Surv(
            Overall.Survival..Months.,
            as.numeric(factor(Overall.Survival.Status)) - 1
        ) ~ group + age_at_index, data=df
    )

}
