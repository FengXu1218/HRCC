## Data preprocessing ------------------------
library(edgeR)
library(pracma)
library(CAFF)
library(Seurat)
library(SeuratDisk)

dir_path <- "./DataSets/nolabel/GSE167607_hFib/"
hFib.data <- Read10X(paste0(dir_path,'rawdata/'))
hFib <- as.data.frame(hFib.data)
# Convert(paste0(dir_path, "filtered_data/velocity_anndata_human_fibroblast_DeepCycle_ISMARA.h5ad"),
#         "h5seurat", overwrite = F, assay = "RNA")

sce <- LoadH5Seurat(paste0(dir_path, "filtered_data/velocity_anndata_human_fibroblast_DeepCycle_ISMARA.h5seurat"))

DimPlot(sce, reduction = "pca", cols = sce@meta.data[["leiden"]])

proliferation_cells <- colnames(sce)[which(sce@meta.data[["leiden"]] == 0)]
proliferation_cells <- gsub("BLKS9:", "",
                            gsub("x", "", proliferation_cells))

hFib <- hFib[, gsub("-1","", colnames(hFib)) %in% proliferation_cells]


# Check the sequencing depth distribution:
quantile(colSums(hFib), probs = c(0,0.1, 0.25,0.5,0.75,0.9, 1))

# Cells with sequencing depth less than 10000 were removed:
colids <- which(colSums(hFib) > 10000)

# Genes with 0 expression in more than 90% of all cells were removed:
rowids <- which(apply(hFib, 1, function(x) length(which(x == 0))/length(x) <= 0.9))

exp <- hFib[rowids, colids]

# cpm transformation：
exp_cpm <- edgeR::cpm(exp, log = T) %>% round(2)
exp_cpm_nolog <- edgeR::cpm(exp, log = F) %>% round(2)

# write.csv(exp_cpm, paste0(dir_path, "cpm/hFib_cpm.csv"))
# write.csv(exp_cpm_nolog, paste0(dir_path, "cpm/hFib_cpm_nolog.csv"))
save(exp, exp_cpm,exp_cpm_nolog, file = paste0(dir_path, "Rdata/exp.Rdata"))

# extracte the expression matrix of cell cycle related genes:
load("Rdata/cc_genes.Rdata")
cc_genes <- intersect(seurat_cc_genes, rownames(exp_cpm))
data_cc <- exp_cpm[cc_genes,]
dim(data_cc)
# write.csv(data_cc, paste0(dir_path, "cpm/hFib_cc_cpm.csv"))

# the pca_plot function to see how the data is distributed on pc1 and pc2:
p_x <- pca_plot(data_cc, time.iqr = 1)
p_x$plot1
p_x$plot2

# remove_pc_outliers
final_index <- p_x$final_index
data_cc_new <- data_cc[, final_index]

write.csv(data_cc_new, paste0(dir_path, "filtered_data/hFib_cc_cpm_filter.csv"))
save(data_cc, data_cc_new, file = paste0(dir_path, "/Rdata/exp_cc_cpm.Rdata"))

