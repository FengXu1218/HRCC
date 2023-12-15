# 1. Quality Control ------------------------------------------------------
library(edgeR)
library(CAFF)
library(pracma)

dir_path <- "./DataSets/label/GSE64016_hESC/"
hESC <- read.csv(paste0(dir_path, "rawdata/GSE64016_H1andFUCCI_normalized_EC.csv"), header = T, row.names = 1)

# Check the sequencing depth distribution:
quantile(colSums(hESC), probs = c(0,0.1, 0.25,0.5,0.75,0.9, 1))

# Cells with sequencing depth less than 10000 were removed:
colids <- which(colSums(hESC) > 10000)

# Genes with 0 expression in more than 70% of all cells were removed:
rowids <- which(apply(hESC, 1, function(x) length(which(x == 0))/length(x) <= 0.7))

exp <- hESC[rowids, colids]

# log transformation：
exp_tpm <- log(exp+1) %>% round(2)
exp_tpm_nolog <- exp %>% round(2)

# write.csv(exp_tpm, paste0(dir_path, "tpm/hESC_tpm.csv"))
# write.csv(exp_tpm_nolog, paste0(dir_path, "tpm/hESC_tpm_nolog.csv"))
save(exp, exp_tpm, exp_tpm_nolog, file = paste0(dir_path, "Rdata/exp.Rdata"))

# extracte the expression matrix of cell cycle related genes:
load("Rdata/cc_genes.Rdata")
cc_genes <- intersect(seurat_cc_genes, rownames(exp_tpm))
data_cc <- exp_tpm[cc_genes,]

# write.csv(data_cc, paste0(dir_path, "tpm/hESC_cc_tpm.csv"))

# the pca_plot function to see how the data is distributed on pc1 and pc2:
p_x <- pca_plot(data_cc, time.iqr = 1)
p_x$plot1
p_x$plot2

# remove_pc_outliers
final_index <- p_x$final_index
data_cc_new <- data_cc[, final_index]

pca_scatter_plot(prcomp(t(data_cc_new)), title = "hESC")
ggsave(paste0(dir_path, "plots/pca_97_Seurat_genes.pdf"), height = 5, width = 5)

write.csv(data_cc_new, paste0(dir_path, "filtered_data/hESC_cc_tpm_filter.csv"))
save(data_cc, data_cc_new, file = paste0(dir_path, "/Rdata/exp_cc_tpm.Rdata"))

