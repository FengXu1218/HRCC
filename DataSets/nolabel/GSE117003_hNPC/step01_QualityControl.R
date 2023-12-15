# 1. Quality Control ------------------------------------------------------
library(edgeR)
library(CAFF)
library(pracma)
library(Seurat)
library(GenomicFeatures)

dir_path <- "./DataSets/nolabel/GSE117003_hNPC/"
hNPC <- Read10X(paste0(dir_path, "rawdata/GSM3267241/"))
hNPC <- as.data.frame(hNPC)

# Check the sequencing depth distribution:
quantile(colSums(hNPC), probs = c(0,0.1, 0.25,0.5,0.75,0.9, 1))

# Cells with sequencing depth less than 10000 were removed:
colids <- which(colSums(hNPC) > 10000)

# Genes with 0 expression in more than 70% of all cells were removed:
rowids <- which(apply(hNPC, 1, function(x) length(which(x == 0))/length(x) <= 0.96))

exp <- hNPC[rowids, colids]

# cpm transformation：
exp_cpm <- edgeR::cpm(exp, log = T) %>% round(2)
exp_cpm_nolog <- edgeR::cpm(exp, log = F) %>% round(2)

# cpm transformation：
# {
#   load(paste0(dir_path,"Rdata/gtf_df.Rdata"))
#   gene_gtf <- gtf_df %>% filter(type == "gene")
#   table(rownames(hNPC) %in% gene_gtf$gene_name)
#   hNPC <- hNPC[rownames(hNPC) %in% gene_gtf$gene_name,]
#   rownames(gene_gtf) <- gene_gtf$gene_name
#   gene_length <- gene_gtf[match(rownames(hNPC),gene_gtf$gene_name),4]
#   #计算
#   kb <- gene_length / 1000
#   rpk <- hNPC / kb
#   cpm <- t(t(rpk)/colSums(rpk) * 1000000)
#   cpm[1:5, 1:5]
#   #log
#   range(cpm)
#   exp_cpm <- log(cpm+1) %>% round(2)
#   exp_cpm_nolog <- cpm %>% round(2)
# }

# write.csv(exp_cpm, paste0(dir_path, "cpm/hNPC_cpm.csv"))
# write.csv(exp_cpm_nolog, paste0(dir_path, "cpm/hNPC_cpm_nolog.csv"))
save(exp, exp_cpm, exp_cpm_nolog, file = paste0(dir_path, "Rdata/exp.Rdata"))

# extracte the expression matrix of cell cycle related genes:
load("Rdata/cc_genes.Rdata")

cc_genes <- intersect(seurat_cc_genes, rownames(exp_cpm))
data_cc <- exp_cpm[cc_genes,]

# write.csv(data_cc, paste0(dir_path, "cpm/hNPC_cc_cpm.csv"))

# the pca_plot function to see how the data is distributed on pc1 and pc2:
p_x <- pca_plot(data_cc, time.iqr = 1)
p_x$plot1
p_x$plot2

# remove_pc_outliers
final_index <- p_x$final_index
data_cc_new <- data_cc[, final_index]
# pdf(paste0(dir_path, "plots/PCA.pdf"), height = 8, width = 10)
# pca_scatter_plot(prcomp(t(data_cc_new)), title = "hNPC")
# dev.off()

write.csv(data_cc_new, paste0(dir_path, "filtered_data/hNPC_cc_cpm_filter.csv"))
save(data_cc, data_cc_new, file = paste0(dir_path, "/Rdata/exp_cc_cpm.Rdata"))

