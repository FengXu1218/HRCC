# 1. Quality Control ------------------------------------------------------
library(edgeR)
library(CAFF)
library(pracma)
library(Seurat)
library(SeuratDisk)

dir_path <- "GSE227555/"
do_QC <- function(dir_name){
  data <- Read10X(dir_name)
  data <- data[["Gene Expression"]]
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
  exp_cpm <- log(cpm(exp)+1) %>% round(2)

  return(exp_cpm)
}

data_UT_cpm <- do_QC(paste0(dir_path, "rawdata/UT"))
data_50_cpm <- do_QC(paste0(dir_path, "rawdata/50"))
# data_75_cpm <- do_QC(paste0(dir_path, "rawdata/75"))

dim(data_UT_cpm)
dim(data_50_cpm)
# dim(data_75_cpm)

save(list = ls(pattern = "cpm$"),
     file = paste0(dir_path, "Rdata/exp_cpm.Rdata"))

############# 提取细胞周期基因 -----------------
dir_path <- "GSE227555/"
load("../PooledAnalysis/Rdata/gene_list.Rdata")
load(paste0(dir_path, "Rdata/exp_cpm.Rdata"))

UT_cpm_cc <- data_UT_cpm[intersect(gene_8, rownames(data_UT_cpm)), ]
dim(UT_cpm_cc)

H50_cpm_cc <- data_50_cpm[intersect(gene_8, rownames(data_50_cpm)), ]
dim(H50_cpm_cc)

# H75_cpm_cc <- data_75_cpm[intersect(gene_8, rownames(data_75_cpm)), ]
# dim(H75_cpm_cc)
