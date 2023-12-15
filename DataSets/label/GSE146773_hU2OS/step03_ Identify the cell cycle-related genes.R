## Find cell cycle-related genes --------------------------------
library(CAFF)

dir_path <- "./DataSets/label/GSE146773_hU2OS/"
load(paste0(dir_path,"/Rdata/exp.Rdata"))
load(paste0(dir_path,"/Rdata/exp_cc_tpm.Rdata"))
load(paste0(dir_path, "/Rdata/pseudo_order.Rdata"))
load("Rdata/cc_genes.Rdata")

# res_cc_genes <- find_cc_genes(exp_tpm[,colnames(data_cc_new)],
#                               pseudo_order, dcor = F,
#                               lag = 1)
# dim(res_cc_genes)

# 300 cells were randomly selected：
tmp <- exp_tpm[,colnames(data_cc_new)][, pseudo_order][, sort(sample(1:ncol(data_cc_new), 300))]
res_cc_genes <- find_cc_genes(tmp, 1:300, dcor = F)
dim(res_cc_genes)

save(res_cc_genes,file = paste0(dir_path, "/Rdata/res_cc_genes.Rdata"))

# cc_genes_tmp <- Reduce(intersect, list(res_cc_genes$cc_gene_test[order(res_cc_genes$p.adj)[1:200]],
#                                        res_cc_genes$cc_gene_test[order(res_cc_genes$dcor_values, decreasing = T)[1:200]],
#                                        res_cc_genes$cc_gene_test[order(res_cc_genes$knn_mi_values, decreasing = T)[1:200]]))

Reduce(intersect, list(res_cc_genes$cc_gene_test,
                       cyclebase3.0_genes,
                       seurat_cc_genes,
                       GO0007049_genes))

cc_gene_list <- list(res_cc_genes$cc_gene_test, seurat_cc_genes,
                     GO0007049_genes, cyclebase3.0_genes)
names(cc_gene_list) <- c("find in hU2OS" , "Seurat" , "GO0007049", "cyclebase3.0")
save_path <- paste0(dir_path, "plots/")

library(VennDiagram)
venn_upset_plot(cc_gene_list, save_path)
