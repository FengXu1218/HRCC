# 1. Quality Control ------------------------------------------------------
library(edgeR)
library(CAFF)
library(pracma)

dir_path <- "./DataSets/label/GSE146773_hU2OS/"
hU2OS <- read.csv(paste0(dir_path, "rawdata/GSE146773_Tpms.csv"), header = T, row.names = 1)
hU2OS <- as.data.frame(t(hU2OS))

# Gene identifier transform
hu_gene_id=rownames(hU2OS)

library(biomaRt)
listEnsembl()
ensembl <- useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl")
datasets <- listDatasets(ensembl)
attributes = listAttributes(ensembl)
hu_gene_name <- getBM(attributes=c("external_gene_name",'ensembl_gene_id'),
                      filters = "ensembl_gene_id",
                      values = hu_gene_id, mart = ensembl)

# remove rows containing NA or duplicates
hU2OS$gene=hu_gene_name[match(rownames(hU2OS),hu_gene_name$ensembl_gene_id),1]
hU2OS <- hU2OS[!is.na(hU2OS$gene),]
hU2OS <- hU2OS[!duplicated(hU2OS$gene),]
rownames(hU2OS) <- hU2OS$gene
hU2OS <- hU2OS[,-ncol(hU2OS)]

# Check the sequencing depth distribution:
quantile(colSums(hU2OS), probs = c(0,0.1, 0.25,0.5,0.75,0.9, 1))

# Cells with sequencing depth less than 10000 were removed:
colids <- which(colSums(hU2OS) > 10000)

# Genes with 0 expression in more than 70% of all cells were removed:
rowids <- which(apply(hU2OS, 1, function(x) sum(x == 0)/length(x) <= 0.7))

exp <- hU2OS[rowids, colids]

exp_tpm <- log2(exp+1)
save(exp, exp_tpm, file = paste0(dir_path, "Rdata/exp.Rdata"))

# extracte the expression matrix of cell cycle related genes:
load("Rdata/cc_genes.Rdata")
cc_genes <- intersect(seurat_cc_genes, rownames(exp_tpm))
data_cc <- exp_tpm[cc_genes,]

# write.csv(data_cc, paste0(dir_path, "tpm/hU2OS_cc_tpm.csv"))

# the pca_plot function to see how the data is distributed on pc1 and pc2:
p_x <- pca_plot(data_cc, time.iqr = 1.25)
p_x$plot1
p_x$plot2

# remove_pc_outliers
final_index <- p_x$final_index
data_cc_new <- data_cc[, final_index]

pca_scatter_plot(prcomp(t(data_cc_new)), title = "hU2OS")

ggsave(paste0(dir_path, "plots/pca_97_Seurat_genes.pdf"), height = 5, width = 7)
write.csv(data_cc_new, paste0(dir_path, "filtered_data/hU2OS_cc_tpm_filter.csv"))
save(data_cc, data_cc_new, file = paste0(dir_path, "/Rdata/exp_cc_tpm.Rdata"))
