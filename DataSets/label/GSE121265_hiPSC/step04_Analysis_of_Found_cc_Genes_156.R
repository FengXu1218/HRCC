# 5. Second round of analysis using new cell cycle-related genes ----------------------------------------------
# 5.1 PCA of cc genes -------------------------------------------------
library(CAFF)
library(cowplot)

dir_path <- "./DataSets/label/GSE121265_hiPSC/"

load(paste0(dir_path, "Rdata/exp.Rdata"))
load("Rdata/gene_list.Rdata")
load("Rdata/cc_genes.Rdata")
# load("Rdata/caff_cc_genes.Rdata")

# pca scatter plot:
# p1 <- pca_scatter_plot(prcomp(t(exp_cpm[intersect(rownames(exp_cpm),
#                                                   seurat_cc_genes),])),
#                        add_legend = F,
#                        title = "hiPSC", sub_title = "97 cell cycle-related genes from Seurat")
# #p1
# p2 <- pca_scatter_plot(prcomp(t(exp_cpm[intersect(rownames(exp_cpm),
#                                                   caff_cc_genes),])),
#                        add_legend = F,
#                        title = "hiPSC",
#                        sub_title = paste0(length(caff_cc_genes),
#                                           " cell cycle-related genes from CAFF"))
#
# p_list <- list(p1,p2)
# n = 3
# for (i in 5:14) {
#   p_list[[n]] <- pca_scatter_plot(prcomp(t(exp_cpm[intersect(rownames(exp_cpm),
#                                                              get(paste0("gene_", i))),])),
#                                   add_legend = F,
#                                   title = "hiPSC", sub_title = paste0(length(get(paste0("gene_", i))),
#                                                                       " cell cycle-related genes among ", i, " times"))
#   n = n+1
# }
# plot_grid(plotlist = p_list, ncol = 4)
# ggsave(paste0(dir_path, "plots/pca_cc_genes.pdf"), height = 11, width = 15)

# extracte the expression matrix of cell cycle related genes:
cc_genes <- intersect(gene_8, rownames(exp_cpm))
data_cc <- exp_cpm[cc_genes,]

# the pca_plot function to see how the data is distributed on pc1 and pc2:
p_x <- pca_plot(data_cc, time.iqr = 1)
p_x$plot1
p_x$plot2

# remove_pc_outliers
final_index <- p_x$final_index
data_cc_new <- data_cc[, final_index]

#usethis::use_data(data_cc, data_cc_new)

## build_pseudo_order --------------------------------
# data(data_cc_new)
## Our method ==============
pseudo_order_list <- build_pseudo_order(data_cc_new, method = "Default")
pseudo_order <- pseudo_order_list$pseudo_order
pseudo_order_rank <- pseudo_order_list$pseudo_order_rank

## Comparison above the three groups ==============
load(paste0(dir_path, "Rdata/pseudo_rank.Rdata"))

rank_data$pseudo_order_rank2 <- pseudo_order_rank[rownames(rank_data)]
rank_data <- rank_data[order(rank_data$pseudo_order_rank),]
angle_rank_data <- rank_data

for (i in 1:ncol(rank_data)) {
  angle_rank_data[,i] <- rank_to_angle(rank_data[,i], max_rank = max(rank_data[,i]))
}

rank_data <- na.omit(rank_data)
p_list <- circ_cor_scatter_plot(rank_data,
                                cor_method = "Spearman_cor")

plot_grid(plotlist = p_list$plot, ncol = 3)
ggsave(paste0(dir_path, "plots/cor_of_pseudo_orders_152.pdf"), height = 8, width = 12)


## Fourier transform filtering ---------------------
## Manually replicate three periodic data for analysis ========
data_cc_ordered <- data_cc_new[, pseudo_order]
save(data_cc_new, data_cc_ordered, pseudo_order,
     file = paste0(dir_path, "Rdata/data_new.Rdata"))

data_cc_bind <- cbind(data_cc_ordered, data_cc_ordered, data_cc_ordered)
X_fftf <- data_cc_bind

for (i in 1:nrow(data_cc_bind)) {
  res <- fftf(1:ncol(data_cc_bind), exp(as.numeric(data_cc_bind[i,])), 1/(ncol(data_cc_bind)/6), plot = F)
  X_fftf[i,] <- res[[1]]
  print(i)
}

## Judge the clockwise direction：
X_fftf_hiPSC <- X_fftf

load("DataSets/nolabel/GSE129447_HeLa1/Rdata/X_fftf2_152.Rdata")
X_fftf_HeLa1 <- X_fftf

X_fftf_names <- c("HeLa1", "hiPSC")
inter_genes <- intersect(rownames(X_fftf_HeLa1), rownames(X_fftf_hiPSC))
angle_data <- as.data.frame(matrix(NA, nrow = length(inter_genes),
                                   ncol = length(X_fftf_names)))
rank_data <- as.data.frame(matrix(NA, nrow = length(inter_genes),
                                   ncol = length(X_fftf_names)))

for (i in 1:length(X_fftf_names)) {
  assign(paste0("rank_", X_fftf_names[i]),
         rank(get_row_order(get(paste0("X_fftf_", X_fftf_names[i]))[inter_genes,])[[1]]))
  assign(paste0("angle_", X_fftf_names[i]),
         rank_to_angle(get(paste0("rank_", X_fftf_names[i])), length(inter_genes)))
  angle_data[,i] <- get(paste0("angle_", X_fftf_names[i]))
  rank_data[,i] <- get(paste0("rank_", X_fftf_names[i]))
}

colnames(angle_data) <- X_fftf_names
colnames(rank_data) <- X_fftf_names
p <- circ_cor_scatter_plot(rank_data, cor_method = "Spearman_cor")
p[["plot"]]

# X_fftf_hiPSC <- X_fftf_hiPSC[,ncol(X_fftf_hiPSC):1]
#
# for (i in 1:length(X_fftf_names)) {
#   assign(paste0("rank_", X_fftf_names[i]),
#          rank(get_row_order(get(paste0("X_fftf_", X_fftf_names[i]))[inter_genes,])[[1]]))
#   assign(paste0("angle_", X_fftf_names[i]),
#          rank_to_angle(get(paste0("rank_", X_fftf_names[i])), length(inter_genes)))
#   angle_data[,i] <- get(paste0("angle_", X_fftf_names[i]))
# }
#
# colnames(angle_data) <- X_fftf_names
# p <- circ_cor_scatter_plot(angle_data, cor_method = "Spearman_cor")
# p[["plot"]]

X_fftf <- X_fftf_hiPSC

# write.csv(X_fftf, paste0(dir_path, "X_FFTF_data/hiPSC_X_fftf2_152.csv"))
write.table(round(t(X_fftf[,(ncol(X_fftf)/3+1):(ncol(X_fftf)*2/3)]),3),
            paste0(dir_path, "X_FFTF_data/hiPSC_X_fftf_152.txt"), sep = ",",
            col.names = F, row.names = F)
save(X_fftf, file = paste0(dir_path, "Rdata/X_fftf2_152.Rdata"))
# usethis::use_data(X_fftf)

## Visualization =============
load(paste0(dir_path, "Rdata/X_fftf2_152.Rdata"))

# Reorder the genes and make a heatmap:
save_path <- paste0(dir_path, "plots/")
p_ind <- heatmap_plot(X_fftf, name = "hiPSC_152_cc_genes", output = T, save_path = save_path)

pca_scatter_plot(prcomp(t(X_fftf[,(ncol(X_fftf)/3+1):(ncol(X_fftf)*2/3)])), title = "X_fftf_hiPSC")
ggsave(paste0(save_path, "pca_fftf_152_genes.pdf"), height = 5, width = 6)

library(ComplexHeatmap)

pdf(paste0(save_path, "Heatmap_fftf_152_genes.pdf"), height = 6, width = 6)
Heatmap(X_fftf[p_ind$final_index,(ncol(X_fftf)/3+1):(ncol(X_fftf)*2/3)],
        cluster_rows = F, cluster_columns = F,
        show_column_names = F, row_names_gp = gpar(fontsize = 3))
dev.off()











