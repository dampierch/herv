library(readr)
library(xtable)

data_dir <- "/scratch/chd5n/herv/annotations/"
targets <- list()
targets[[1]] <- paste0(data_dir, "pheno_A.tsv")
targets[[2]] <- paste0(data_dir, "pheno_C.tsv")
targets[[3]] <- paste0(data_dir, "pheno_B.tsv")

dfs <- list()
for (i in seq_len(3)) {
    dfs[[i]] <- readr::read_tsv(targets[[i]])
}

lapply(dfs, fun <- function(df) {
    summary(factor(df$race))
})
