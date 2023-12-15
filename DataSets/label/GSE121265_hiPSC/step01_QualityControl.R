# 1. Quality Control ------------------------------------------------------
library(edgeR)
library(org.Hs.eg.db)
library(CAFF)
library(pracma)
library(clusterProfiler)

dir_path <- "./DataSets/label/GSE121265_hiPSC/"

data <- readRDS(paste0(dir_path, "rawdata/GSE121265_eset-final.rds"))

iPSC <- as.data.frame(data@assayData[["exprs"]][-1,])
ids <- bitr(rownames(iPSC), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
ids <- ids[!duplicated(ids$ENSEMBL), ]
ids <- ids[!duplicated(ids$SYMBOL), ]
iPSC$ENSEMBL <- rownames(iPSC)
iPSC_new <- inner_join(iPSC, ids)
rownames(iPSC_new) <- iPSC_new$SYMBOL
iPSC_new <- iPSC_new[,-c(889:890)]
iPSC <- iPSC_new

# Check the sequencing depth distribution:
quantile(colSums(iPSC), probs = c(0,0.1, 0.25,0.5,0.75,0.9, 1))

# Cells with sequencing depth less than 10000 were removed:
colids <- which(colSums(iPSC) > 10000)

# Genes with 0 expression in more than 70% of all cells were removed:
rowids <- which(apply(iPSC, 1, function(x) length(which(x == 0))/length(x) <= 0.7))

exp <- iPSC[rowids, colids]

# CPM:
exp_cpm <- edgeR::cpm(exp, log = T)
exp_cpm_nolog <- edgeR::cpm(exp, log = F)

# write.csv(exp_cpm, paste0(dir_path, "cpm/iPSC_cpm.csv"))
# write.csv(exp_cpm_nolog, paste0(dir_path, "cpm/iPSC_cpm_nolog.csv"))
save(exp, exp_cpm, exp_cpm_nolog, file = paste0(dir_path, "Rdata/exp.Rdata"))

# extracte the expression matrix of cell cycle related genes:
load("Rdata/cc_genes.Rdata")
cc_genes <- intersect(seurat_cc_genes, rownames(exp_cpm))
data_cc <- exp_cpm[cc_genes,]

# write.csv(data_cc, paste0(dir_path, "cpm/iPSC_cc_cpm.csv"))

# the pca_plot function to see how the data is distributed on pc1 and pc2:
p_x <- pca_plot(data_cc, time.iqr = 1)
p_x$plot1
p_x$plot2

# remove_pc_outliers
final_index <- p_x$final_index
data_cc_new <- data_cc[, final_index]

pca_scatter_plot(prcomp(t(data_cc_new)), title = "iPSC")

ggsave(paste0(dir_path, "plots/pca_97_Seurat_genes.pdf"), height = 5, width = 5)
write.csv(data_cc_new, paste0(dir_path, "filtered_data/iPSC_cc_cpm_filter.csv"))
save(data_cc, data_cc_new, file = paste0(dir_path, "/Rdata/exp_cc_cpm.Rdata"))
