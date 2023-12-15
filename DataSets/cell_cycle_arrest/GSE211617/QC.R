# 1. Quality Control ------------------------------------------------------
library(edgeR)
library(CAFF)
library(pracma)
library(Seurat)
library(SeuratDisk)

dir_path <- "GSE211617/"

do_QC <- function(dir_name){
  data <- Read10X(dir_name)

  data_sce <- CreateSeuratObject(data, min.cells = 3, min.features = 200)
  dim(data_sce)

  # Check the sequencing depth distribution:
  quantile(colSums(data_sce), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))

  # Cells with sequencing depth less than 10000 were removed:
  colids <- which(colSums(data_sce) > 10000)

  # Genes with 0 expression in more than 70% of all cells were removed:
  # rowids <- which(apply(data_sce, 1, function(x) length(which(x == 0))/length(x) <= 0.9))

  data_sce2 <- data_sce[, colids]
  exp <- as.matrix(data_sce2@assays[["RNA"]]@counts)
  dim(exp)

  # cpm transformation：
  exp_cpm <- cpm(exp, log = T) %>% round(2)

  return(exp_cpm)
}

data_IR2h_cpm <- do_QC(paste0(dir_path, "rawdata/IR_2h"))
data_IR6h_cpm <- do_QC(paste0(dir_path, "rawdata/IR_6h"))
data_Ctrl_cpm <- do_QC(paste0(dir_path, "rawdata/Ctrl_P"))

dim(data_IR2h_cpm)
dim(data_IR6h_cpm)
dim(data_Ctrl_cpm)

save(list = ls(pattern = "cpm$"),
     file = paste0(dir_path, "Rdata/exp_cpm.Rdata"))
