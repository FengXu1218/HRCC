## build_pseudo_order --------------------------------
library(CAFF)

dir_path <- "./DataSets/nolabel/GSE167607_hFib/"
load(paste0(dir_path, "/Rdata/exp_cc_cpm.Rdata"))

## Our method ==============
pseudo_order_list <- build_pseudo_order(data_cc_new, method = "Default", reverse = F)
pseudo_order <- pseudo_order_list$pseudo_order
pseudo_order_rank <- pseudo_order_list$pseudo_order_rank

## tricycle method ==============
tricycle_order_rank <- build_pseudo_order(data_cc_new, method = "tricycle")


# ## reCAT method ==============
# data_cc_cyclebase <- exp_cpm[intersect(cyclebase3.0_genes, rownames(exp_cpm)),
#                              intersect(colnames(data_cc_new), colnames(exp_cpm))]
# dim(data_cc_cyclebase)
#
# reCAT_order <- build_pseudo_order(data_cc_cyclebase, method = "reCAT")
# reCAT_order_rank <- reCAT_order[[2]]
#
## Comparison above the three groups ==============
rank_data <- data.frame(pseudo_order_rank = pseudo_order_rank,
                        tricycle_order_rank = tricycle_order_rank #,
                        # reCAT_order_rank = reCAT_order_rank
)
rank_data <- rank_data[order(rank_data$pseudo_order_rank),]

p_list <- circ_cor_scatter_plot(rank_data, cor_method = "Spearman_cor")

library(cowplot)
p <- plot_grid(plotlist = p_list$plot, ncol = 1)
p

ggsave(paste0(dir_path, "/plots/cor_of_pesudo_orders.pdf"),
       plot = p, height = 4, width = 4)
save(rank_data, file = paste0(dir_path, "/Rdata/pseudo_rank.Rdata"))
save(pseudo_order, file = paste0(dir_path, "/Rdata/pseudo_order.Rdata"))



