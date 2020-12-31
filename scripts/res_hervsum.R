## to be sourced in results.R


count_herv_ids <- function(gtf, ids) {
    cat("Counting total herv transcripts\n")
    cmd = paste(
        paste("zcat", gtf),
        "awk -v FS=\"\\t\" '!($0~/^#/) && $3 == \"transcript\" {print $0}'",
        "wc -l",
        sep=" | "
    )
    x1 <- as.numeric(system(cmd, intern=TRUE))
    cat("Counting select herv ids\n")
    df <- data.frame(readr::read_tsv(ids))
    x2 <- nrow(df)
    x3 <- length(unique(df$herv_id))
    nms <- c("Total Transcripts", "Select Transcripts", "HERV Loci")
    return(setNames(list(x1, x2, x3), nms))
}


count_herv_cats <- function(cv) {
    nms <- base::unlist(
        base::lapply(
            base::lapply(
                strsplit(cv, "[|]"), unpack_id
            ),
            head, 1
        ),
        use.names=FALSE
    )
    cats <- c(
        "HERVE", "HERVH", "HERVK", "HERVL"
    )
    cnts <- list()
    for (e in cats) {
        cntr <- 0
        for (i in nms) {
            if (i == e) {
                cntr <- cntr + 1
            }
        }
        cnts[[e]] <- cntr
    }
    cnts[["Other"]] <- length(nms) - sum(unlist(cnts))
    return(cnts)
}


fill_params <- function(l, params) {
    if ("Total Transcripts" %in% names(l)) {
        params[["ggp_title"]] <- "Number of HERV Elements Analyzed"
        n <- format(l[["Total Transcripts"]], big.mark=",")
        params[["ggp_sub"]] <- paste("From", n, "starting transcripts")
        l <- l[names(l) != "Total Transcripts"]
    } else {
        params[["ggp_title"]] <- "Groups of HERV Loci Tested"
        params[["ggp_sub"]] <- element_blank()
    }
    params[["ggp_ylab"]] <- "Counts"
    return(setNames(list(l, params), c("l", "params")))
}


format_data <- function(l) {
    df <- data.frame(base::t(data.frame(l, check.names=FALSE)))
    df[, 2] <- rownames(df)
    colnames(df) <- c("counts", "categories")
    if ("Other" %in% rownames(df)) {
        tmp <- df[df$categories != "Other", ]
        df$categories <- factor(
            df$categories,
            levels=c(
                rownames(tmp[order(tmp$counts, decreasing=TRUE), ]),
                "Other"
            )
        )
    } else {
        df$categories <- factor(
            df$categories,
            levels=rownames(df[order(df$counts, decreasing=TRUE), ])
        )
    }
    return(df)
}


plot_hrvsum <- function(inlist) {
    ggp_l <- list()
    cntr <- 1
    for (l in inlist) {
        params <- list()
        int_l <- fill_params(l, params)
        df <- format_data(int_l$l)
        p <- int_l$params
        ggp_l[[cntr]] <- ggplot(df, aes(x=categories, y=counts)) +
            geom_col() +
            labs(title=p[["ggp_title"]], subtitle=p[["ggp_sub"]],
                x=element_blank(), y=p[["ggp_ylab"]]
            ) +
            geom_text(
                stat='identity',
                aes(label=format(counts, big.mark=",")),
                vjust=-0.5,
                size=3
            ) +
            scale_y_continuous(limits=c(0, max(df$counts) * 1.05)) +
            ggp_theme_col
        cntr <- cntr + 1
    }
    return(ggp_l)
}


write_hrvsum <- function(plotlist, target) {
    p_width <- 7
    p_height <- 3.5
    pdf(file=target, width=p_width, height=p_height)
    print(
        cowplot::plot_grid(
            plotlist=plotlist,
            nrow=1,
            ncol=2,
            labels=c("A", "B")
        )
    )
    dev.off()
    cat("figure written to", target, "\n")
}
