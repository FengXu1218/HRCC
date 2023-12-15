# 3.Find cell cycle-related genes --------------------------------
library(CAFF)

dir_path <- "./DataSets/nolabel/GSE117003_hNPC/"
load(paste0(dir_path,"/Rdata/exp.Rdata"))
load(paste0(dir_path,"/Rdata/exp_cc_cpm.Rdata"))
load(paste0(dir_path, "/Rdata/pseudo_order.Rdata"))
load("Rdata/cc_genes.Rdata")

# 3.1 Identify the cell cycle-related genes -------------------------------
res_cc_genes <- find_cc_genes(exp_cpm[, colnames(data_cc_new)],
                              pseudo_order, dcor = F)
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
names(cc_gene_list) <- c("find in hNPC" , "Seurat" , "GO0007049", "cyclebase3.0")

library(VennDiagram)
save_path <- paste0(dir_path, "plots/")

venn_upset_plot(cc_gene_list, save_path)
