###################### Figure 1全部图形 ####################
library(ggplot2)

####### Fig 1A: PCA of Seurat in Hela1 -----------
dir_path <- "./DataSets/nolabel/GSE129447_HeLa1/"

source("R/visualization.R")
source("R/build_pseudo_order.R")
source("R/data_preprocessing.R")

load(paste0(dir_path, "Rdata/exp.Rdata"))
load("Rdata/cc_genes.Rdata")
cc_genes <- intersect(seurat_cc_genes, rownames(exp_cpm))
data_cc <- exp_cpm[cc_genes,]

# the pca_plot function to see how the data is distributed on pc1 and pc2:
p_x <- pca_plot(data_cc, time.iqr = 1)
p_x$plot1

# remove_pc_outliers
final_index <- p_x$final_index
data_cc_new <- data_cc[, final_index]

library(RColorBrewer)
library(pracma)
p_A <- pca_scatter_plot(prcomp(t(data_cc_new)), title = "PCA of Seurat CC Genes in Hela1",
                       # add_legend = F,
                       # colors = colorRampPalette(c("#2d3d98","#46a6dd","#a4d283",
                       #                           "#feeb32", "#b01f24","#2d3d98"))(10)
)

####### Fig 1B: cor of pseudo time ---------
load("DataSets/nolabel/GSE129447_HeLa1/Rdata/pseudo_rank.Rdata")

p_list <- cor_scatter_plot(rank_data)

library(cowplot)
p_B <- plot_grid(plotlist = p_list, ncol = 3)


######### Fig 1C：bar plot =============
load("Rdata/gene_list.Rdata")

gene_names <- ls(pattern = "gene_[0-9].*")
len_gene_list <- c()

for (i in 1:length(gene_names)) {
  len_gene_list[i] <- length(get(gene_names[i]))
}

len_gene_list <- sort(len_gene_list, decreasing = T)

bar_data <- data.frame("Times_of_Repetition" = 4:12,
                       "Number_of_Genes" = len_gene_list)

p_C <- ggplot(bar_data) +
  geom_col(aes(Times_of_Repetition, Number_of_Genes,
               fill = Times_of_Repetition)) +
  geom_text(aes(Times_of_Repetition, Number_of_Genes+15,
                label = Number_of_Genes))+
  scale_fill_gradientn(colors = c("#3d4ba0","#bb6bab", "#e1cf5b",
                                  "#499a44", "#3d4ba0")) +
  scale_x_continuous(breaks = 4:14)+
  xlab("Times of Repetition") +
  ylab("Number of Genes") +
  theme_classic() +
  theme(legend.position = "none")

############ fig 1d: Venn Diagram ------------
library(VennDiagram)

cc_gene_list <- list(gene_8, seurat_cc_genes,
                     GO0007049_genes, cyclebase3.0_genes)
names(cc_gene_list) <- c("Our featured genes" , "Seurat" ,
                         "GO0007049", "Cyclebase3.0")
category.names <- names(cc_gene_list)
fill_color <- brewer.pal(length(cc_gene_list), "Set3")
save_path = "./Figures/figure1/"
venn_name = "fig_1d"

venn.diagram(x = cc_gene_list,
             disable.logging = T,
             category.names = category.names,
             filename = paste0(save_path, venn_name, ".tiff"),
             imagetype = "tiff",
             height = 1000, width = 1000, resolution = 300, compression = "lzw",
             col = "white", lty = 1, lwd = 1, fill = fill_color,
             cat.pos = c(-5, 0, 0, 5),
             alpha = 0.9, label.col = "black", cex = 0.5, fontfamily = "serif",
             fontface = "bold", cat.col = "black", cat.cex = 0.6,
             cat.fontfamily = "serif")







