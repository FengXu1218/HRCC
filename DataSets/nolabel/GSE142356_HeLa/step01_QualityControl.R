# 1. Quality Control ------------------------------------------------------
library(edgeR)
library(CAFF)
library(pracma)

dir_path <- "./DataSets/nolabel/GSE142356_HeLa/"
HeLa <- read.table(paste0(dir_path, "rawdata/GSM4226257_out_gene_exon_tagged.dge_exonsds_046.txt"), header = T, row.names = 1)

# Check the sequencing depth distribution:
quantile(colSums(HeLa), probs = c(0,0.1, 0.25,0.5,0.75,0.9, 1))

# Cells with sequencing depth less than 10000 were removed:
colids <- which(colSums(HeLa) > 10000)

# Genes with 0 expression in more than 70% of all cells were removed:
rowids <- which(apply(HeLa, 1, function(x) length(which(x == 0))/length(x) <= 0.95))

exp <- HeLa[rowids, colids]

# cpm transformation：
exp_cpm <- edgeR::cpm(exp, log = T) %>% round(2)
exp_cpm_nolog <- cpm(exp, log = F) %>% round(2)

# write.csv(exp_cpm, paste0(dir_path, "cpm/HeLa_cpm.csv"))
# write.csv(exp_cpm_nolog, paste0(dir_path, "cpm/HeLa_cpm_nolog.csv"))
save(exp, exp_cpm, exp_cpm_nolog, file = paste0(dir_path, "Rdata/exp.Rdata"))

# extracte the expression matrix of cell cycle related genes:
load("Rdata/cc_genes.Rdata")
cc_genes <- intersect(seurat_cc_genes, rownames(exp_cpm))
data_cc <- exp_cpm[cc_genes,]

# write.csv(data_cc, paste0(dir_path, "cpm/HeLa_cc_cpm.csv"))

# the pca_plot function to see how the data is distributed on pc1 and pc2:
p_x <- pca_plot(data_cc, time.iqr = 1)
p_x$plot1
p_x$plot2

# remove_pc_outliers
final_index <- p_x$final_index
data_cc_new <- data_cc[, final_index]

# write.csv(data_cc_new, paste0(dir_path, "filtered_data/HeLa_cc_cpm_filter.csv"))
save(data_cc, data_cc_new, file = paste0(dir_path, "/Rdata/exp_cc_cpm.Rdata"))

