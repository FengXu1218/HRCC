####### Figure S2 ---------
load("DataSets/nolabel/GSE129447_HeLa1/Rdata/pseudo_rank.Rdata")

p_list <- cor_scatter_plot(rank_data)

library(cowplot)
p_B <- plot_grid(plotlist = p_list, ncol = 3)
