library(propr)
library(compositions)
library(tidyverse)
library(foreach)
library(doParallel)
cl <- makeCluster(1)
registerDoParallel(cl)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args[1] <- "."
}

#### input counts table with dummpy genes
RIBO <- read.csv(paste(args[1], "/ribo_paired_count_dummy.csv", sep = ""), row.names = 1)#[1:1000, ]
RNA <- read.csv(paste(args[1], "/rna_paired_count_dummy.csv", sep = ""), row.names = 1)#[1:1000, ]
RIBO <- t(RIBO)
RNA <- t(RNA)

print(dim(RIBO))

TE_clr <- function(RIBO, RNA) {
  ### data processing, transfer count to clr
  ### clr normalization for ribo-seq and RNA-seq
  pr_RIBO <- propr(RIBO, metric = "rho", ivar = "clr", alpha = NA, p = 100)
  pr_RNA <- propr(RNA, metric = "rho", ivar = "clr", alpha = NA, p = 100)
  print("transfer data from clr to ilr")
  ### transfer data from clr to ilr
  ### This transformation is crucial as it allows the compositional data 
  ### to be decomposed into an array of uncorrelated variables 
  ### while preserving relative proportions. 
  RIBO_ilr <- clr2ilr(pr_RIBO@logratio)
  RNA_ilr <- clr2ilr(pr_RNA@logratio)
  RIBO_ilr <- as.data.frame(t(RIBO_ilr))
  RNA_ilr <- as.data.frame(t(RNA_ilr))
  # out <- data.frame(matrix(ncol = 0, nrow = nrow(RIBO_ilr) + 1))
  ### calculate proportional regression
  # for (i in 1:ncol(RIBO_ilr)) {
  #   print(i)
  #   m <- summary(lm(RIBO_ilr[, i] ~ RNA_ilr[, i]))
  #   out[, i] <- data.frame(as.numeric(ilr2clr(resid(m))))
  # }
  print("calculate proportional regression")
  out <- foreach(i = 1:ncol(RIBO_ilr), .combine = "cbind", .packages = c("compositions")) %dopar% {
    ### compositional linear regression
    m <- summary(lm(RIBO_ilr[, i] ~ RNA_ilr[, i]))
    ### define the residuals as TE and transfer the data back to clr
    data.frame(as.numeric(ilr2clr(resid(m))))
  }

  colnames(out) <- rownames(RIBO)
  rownames(out) <- colnames(RIBO)
  return(out)
}
human_TE <- TE_clr(RIBO, RNA)
save(human_TE, file = paste(args[1], "/human_TE_sample_level.rda", sep = ""))

infor <- read.csv("data/infor_filter.csv")
### merge the TE based on the cell lines and tissues
df <- merge(infor, t(human_TE), by.x = "experiment_alias", by.y = 0)
df_cell_line <- df %>%
  group_by(cell_line) %>%
  summarize(across(where(is.numeric), mean))
df_cell_line_fin <- data.frame(df_cell_line[, 3:ncol(df_cell_line)])
rownames(df_cell_line_fin) <- df_cell_line$cell_line
write.csv(df_cell_line_fin, paste(args[1], "/human_TE_cellline_all.csv", sep = ""))
stopCluster(cl)
