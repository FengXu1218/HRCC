###################### Figure 2全部图形 ####################
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplotify)
library(cowplot)
library(tidyverse)

####### Fig 2A: HeLa1的傅里叶滤波后热图 -----------
dir_path <- "DataSets/nolabel/GSE129447_HeLa1/"

# source("R/visualization.R")
source("R/fftf.R")

load("Rdata/fftf_6.Rdata")
load("Rdata/fftf_all.Rdata")
load(paste0(dir_path, "Rdata/data_new.Rdata"))
load("Rdata/stage_all_correct.Rdata")

col_fun <- colorRamp2(c(-0.1, 0.2, 0.6, 1),
                      c("#74c464", "#badfb7", "#fdedf6","#f064af"))

HeLa1_fftf_6 <- HeLa1_fftf_6[get_row_order(cbind(HeLa1_fftf_6,
                                                 HeLa1_fftf_6,
                                                 HeLa1_fftf_6))[[2]],]

HeLa1_raw <- data_cc_ordered[rownames(HeLa1_fftf_6), colnames(HeLa1_fftf_6)]

# HeLa1_fftf_6 <- HeLa1_fftf_6[,c(238:282,1:237)]
marker_label <- c("CCNE2", "PCNA", "RRM2",
                  "CDK1", "CCNF", "TOP2A", "TPX2")

ha = rowAnnotation(foo = anno_mark(at = which(rownames(HeLa1_fftf_6) %in% marker_label),
                                   labels = marker_label,
                                   labels_gp = gpar(fontsize = unit(5, "pt"),
                                                    fontface = "bold.italic")))

p_A <- Heatmap(HeLa1_fftf_6,
        # 设置颜色：
        col = col_fun,
        show_column_names =F,
        show_row_names = F,
        cluster_rows = F, cluster_columns = F,
        # column_split = cluster_order,
        # # 调整行标签和列标签的大小：
        # row_names_gp = gpar(fontsize = 3,
        #                     fontface = "italic" # 调整字体为斜体
        # ),
        right_annotation = ha,
        # 调整图例：
        # show_heatmap_legend = FALSE,
        # colorbar:
        # top_annotation = cluster_bar,
        column_title = "Pseudo time",
        column_title_side = "bottom",
        column_title_gp = gpar(fontsize = unit(7, "pt"), fontface = "bold.italic")
)
p_A

####### Fig 2B: 基因峰值出现位置的相关性热图 -----------
load("Rdata/max_index_norm.Rdata")

# 矩阵散点图：
library(Hmisc)

apply(gene_rank_data, 2, function(x) length(which(!is.na(x))))

cor_mat <- matrix(NA, 12, 12)
p_mat <- matrix(NA, 12, 12)
rownames(cor_mat) <- colnames(gene_rank_data)
colnames(cor_mat) <- colnames(gene_rank_data)
rownames(p_mat) <- colnames(gene_rank_data)
colnames(p_mat) <- colnames(gene_rank_data)

for (i in 1:ncol(gene_rank_data)) {
  for (j in 1:ncol(gene_rank_data)) {
    tmp_data <- na.omit(gene_rank_data[, c(i,j)])
    res_cor <- rcorr(tmp_data[,1], tmp_data[,2], type = "pearson")
    cor_mat[i,j] <- res_cor[["r"]][1,2]
    p_mat[i,j] <- res_cor[["P"]][1,2]
  }
}

cor_mat <- round(cor_mat,2)
p_mat[p_mat == 0] <- 1e-15
p_mat[upper.tri(cor_mat)] <- cor_mat[upper.tri(cor_mat)]

# 数据转换：
cor_mat_long <- cor_mat %>%
  as.data.frame() %>%
  mutate(x = factor(rownames(cor_mat), levels = rownames(cor_mat))) %>%
  pivot_longer(cols = !x, names_to = "y", values_to = "cor") %>%
  mutate(y = factor(y, levels = rownames(cor_mat)))

head(cor_mat_long)

# 取相关性矩阵上三角：
cor_mat_up <- cor_mat
cor_mat_up[lower.tri(cor_mat_up, diag = T)] <- NA

cor_mat_up_long <- cor_mat_up %>%
  as.data.frame() %>%
  mutate(x = factor(rownames(cor_mat_up), levels = rev(rownames(cor_mat)))) %>%
  pivot_longer(cols = !x, names_to = "y", values_to = "cor") %>%
  mutate(y = factor(y, levels = rownames(cor_mat)))

# 取p值矩阵下三角：
p_mat_lower <- p_mat
p_mat_lower[!lower.tri(p_mat_lower, diag = F)] <- 1

p_mat_lower_long <- p_mat_lower %>%
  as.data.frame() %>%
  mutate(x = factor(rownames(p_mat_lower), levels = rev(rownames(cor_mat)))) %>%
  pivot_longer(cols = !x, names_to = "y", values_to = "p") %>%
  mutate(y = factor(y, levels = rownames(cor_mat)))

# 显著性标记：
p_mat_lower_long$p_sig <- as.character(symnum(p_mat_lower_long$p,
                                              cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                              symbols = c("***", "**", "*", "")))

ggplot() +
  # p值矩阵：
  geom_point(data = p_mat_lower_long,
             aes(x, y, color = -log10(p)),
             size = 10,
             shape = 15)+
  geom_text(data = p_mat_lower_long,
            aes(x, y, label = p_sig))+
  # 配色：
  scale_color_gradientn(colours = c("#ffffff", "#badfb7", "#74c464"))+
  # 气泡图：
  geom_point(data = cor_mat_up_long, aes(x, y, size = cor, fill = cor),
             shape = 21, color = "#c4bcba")+
  geom_text(data = cor_mat_up_long,
            aes(x, y, label = cor),
            color = "black", size = 2)+
  scale_size_continuous(range = c(3, 6))+
  # 配色：
  scale_fill_gradient2(low = "#ffffff", mid = "#e080af", high = "#f064af")+
  # x轴和y轴扩展
  scale_x_discrete(expand = c(0.043, 0.043))+
  scale_y_discrete(expand = c(0.043, 0.043), position = "right")+
  # 框线
  geom_vline(aes(xintercept =seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  geom_hline(aes(yintercept=seq(0.5, 20.5, 1)), color = "#bbbbbb")+
  # 主题：
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, face = "bold.italic"),
        axis.text.y = element_text(face = "bold.italic"),
        legend.position = "top",
        legend.key.width = unit(0.6, "cm"),
        legend.key.height = unit(0.2, "cm"))+
  guides(fill = guide_colorbar("correlation", title.position = "top",
                               title.theme = element_text(face = "bold.italic")),
         color = guide_colorbar("-log10(P value)", title.position = "top",
                                title.theme = element_text(face = "bold.italic")),
         size = "none")

############### Fig 2C：基因方差图 -------------------
library(cowplot)

# 探索方差的情况
var_all_genes <- apply(max_index_norm, 1, function(x) var(na.omit(x)))

gene_var <- data.frame(Genes = rownames(max_index_norm),
                       Variance = var_all_genes)

MultiplePeaks <- names(var_all_genes)[which(var_all_genes < 0.02 & var_all_genes > 0.005)]
BroadPeaks <- names(var_all_genes)[which(var_all_genes < 0.005 & var_all_genes > 0.002)]
StablePeaks <- names(var_all_genes)[which(var_all_genes < 0.002)]

# 绘图：
Peak_names <- ls(pattern = "Peaks")
save(list = Peak_names, file = "Rdata/peak_genes_3Pattern.Rdata")

##### 每组随机抽10个展示
set.seed(4)

for (j in 1:length(Peak_names)) {
  tmp <- max_index_norm[get(Peak_names[j]), -9]
  tmp <- tmp[levels(reorder(get(Peak_names[j]), density_data[which(density_data$gene %in% get(Peak_names[j])),]$x)), ]
  tmp <- tmp[sort(sample(1:nrow(tmp), 10)),]
  p_list <- list()

  p <- ggplot()+
    geom_density(aes(na.omit(tmp[1,])), color = "#f064af",
                 fill = "#f064af", alpha= 0.6, bw = 0.03)+
    ylab(rownames(tmp)[1])+
    # ylab("")+
    xlab("")+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5))+
    xlim(0,1)
  # p

  p_list[[1]] <- p

  color <- colorRampPalette(c("#f064af", "#e4daa3", "#74c464"))(nrow(tmp))

  for (i in 2:nrow(tmp)) {
    tmp2 <- as.data.frame(t(tmp))
    colnames(tmp2)[i] <- "Pseudo Time"
    p_list[[i]] <- ggplot(tmp2)+
      geom_density(aes(x = `Pseudo Time`),
                   color = color[i],bw = 0.03,
                   fill = color[i], alpha= 0.6)+
      ylab(rownames(tmp)[i])+
      # ylab("")+
      xlab("")+
      theme_minimal()+
      theme(panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5))+
      xlim(0,1)
  }

  p_all <- plot_grid(plotlist = p_list, ncol = 1)

  assign(paste0("p", Peak_names[j]), p_all)

  # ggsave(plot = p_all, filename = paste0("plots/p_", Peak_names[j], ".pdf"),
  #        height = nrow(tmp)*2/3, width = 5)
}

p <- plot_grid(pStablePeaks, pBroadPeaks, pMultiplePeaks, ncol = 3)
p

############### Fig 2D 三组基因的富集分析图 -------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

load("Rdata/peak_genes_3Pattern.Rdata")

enrich_GO <- function(gene_list, filename){
  ego <- enrichGO(gene          = gene_list,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)

  write.csv(ego, paste0("Tables/", filename, "ego.csv"), row.names =F)
  ego <- as.data.frame(ego)
  return(ego)
}

ego_Stable <- enrich_GO(StablePeaks, "StablePeaks_")
ego_Broad <- enrich_GO(BroadPeaks, "BroadPeaks_")
ego_Multiple <- enrich_GO(MultiplePeaks, "MultiplePeaks_")

ego_Stable_filter <- ego_Stable %>%
  slice_min(qvalue, n = 10)

ego_Broad_filter <- ego_Broad %>%
  slice_min(qvalue, n = 10)

ego_Multiple_filter <- ego_Multiple %>%
  slice_min(qvalue, n = 10)

ego_filter <- rbind(ego_Stable_filter, ego_Broad_filter, ego_Multiple_filter)
ego_filter$cluster <- factor(rep(c("Stable", "Broad", "Multiple"), each = 10),
                             levels = c("Stable", "Broad", "Multiple"))
ego_filter$Description <- paste(ego_filter$ID, ego_filter$Description)
ego_filter$Description <- factor(ego_filter$Description,
                                 levels = unique(ego_filter$Description[order(ego_filter$cluster, ego_filter$Count, decreasing = T)]))

shorten_names <- function(x, n_word=10, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x)>n_char))
  {
    if (nchar(x) >n_char) x <- substr(x, 1, n_char)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  }
  else
  {
    return(x)
  }
}

ego_filter$Description <- as.factor(ego_filter$Description)
labels=(sapply(
  levels(ego_filter$Description)[as.numeric(ego_filter$Description)],
  shorten_names))


p_D <- ggplot(ego_filter, aes(cluster, Description)) +
  geom_point(aes(fill=-log10(qvalue), size=Count), shape=21)+
  scale_size_continuous(range = c(4, 10))+
  scale_y_discrete(labels=labels) +
  theme_bw()+
  theme(axis.text.x=element_text(angle=270, hjust = 0, vjust = 0.5),
        axis.text = element_text(face = "bold.italic", color = 'black', size = 10),
        legend.position = "top")+
  scale_fill_gradientn(colors = c("#74c464", "#fcf6f9", "#f064af"))+
  labs(x=NULL,y=NULL)+
  coord_flip()

p_D

########### fig2E&F：功能模块方差图 ---------------------
library(readxl)
library(tidyverse)
library(ggsignif)
library(ggplot2)
library(RColorBrewer)
library(rstatix)

FM_data <- read_excel("Tables/function_module.xlsx", sheet = 2)

max_index_norm <- max_index_norm[,-8]

# 基因所属时期数据 -- 来源于Cell文献
gene_meta <- read_excel("Tables/function_module.xlsx", sheet = 3)

S_genes <- intersect(gene_meta$S, rownames(max_index_norm))
G2M_genes <- intersect(gene_meta$G2M, rownames(max_index_norm))

# 每个数据集都能算一个模块基因方差和整个S期基因方差：
var_S <- apply(max_index_norm[S_genes, ], 2, function(x) var(na.omit(x)))
var_G2M <- apply(max_index_norm[G2M_genes, ], 2, function(x) var(na.omit(x)))

## 构造不同通路基因的方差图数据 ---------
var_FM <- matrix(NA, nrow = nrow(FM_data), ncol = ncol(max_index_norm))

for (i in 1:nrow(FM_data)) {
  FM_name <- FM_data$Description[i]
  FM_genes <- str_split(FM_data[i,]$geneID, "/", simplify = T)

  max_index_FM <- max_index_norm[FM_genes, ]
  var_FM[i, ] <- apply(max_index_FM, 2, function(x) var(na.omit(x)))
}

var_all <- as.data.frame(rbind(rbind(var_S, var_G2M), var_FM))
rownames(var_all) <- c("All genes in S phase",
                       "All genes in G2M phase",
                       FM_data$Description)
colnames(var_all) <- gsub("_5", "", colnames(max_index_norm))

var_all$FM_names <- rownames(var_all)
var_all$phase <- c("S", "G2M", FM_data$Phase)
var_all_long <- var_all %>%
  pivot_longer(-c(FM_names, phase)) %>%
  na.omit() %>%
  group_by(FM_names) %>%
  mutate(med = median(value),
         low = quantile(value, 0.25),
         high = quantile(value, 0.75)) %>%
  ungroup()


## 绘图 --------------
var_all_S <- var_all_long %>% filter(phase == "S")

stat_t.test_S <- var_all_S %>%
  t_test(value ~ FM_names, ref.group = "All genes in S phase")

p_S <- ggplot(var_all_S)+
  geom_errorbar(aes(x = reorder(FM_names, med, decreasing = T),
                    color = FM_names, ymin = low, ymax = high), width = 0.2)+
  geom_bar(aes(x = reorder(FM_names, med, decreasing = T),
           fill = FM_names, color = FM_names,
           y = med), data = unique(var_all_S[,c(1,5:7)]),
           width = 0.8, alpha = 0.5,
           stat = "identity")+
  geom_text(aes(reorder(FM_names, med, decreasing = T),
                y = 0.0002,
                label = format(med, scientific=TRUE, digits=2)),
            size = 2)+
  geom_text(data = stat_t.test_S[match(levels(reorder(var_all_S$FM_names, var_all_S$med,
                                                      decreasing = T))[-1],
                                       stat_t.test_S$group2),],
            aes(x = 2:7, y = 0.010, label = p.adj.signif))+
  scale_y_continuous(position = "right")+
  # 颜色模式：
  scale_fill_manual(values = c("#dd67a7", "#f0ad9c", "#85c2a3",
                               "#aec221", "#70be6f", "#ef9c1e", "#17615e"))+
  scale_color_manual(values = c("#dd67a7", "#f0ad9c", "#85c2a3",
                                "#aec221", "#70be6f", "#ef9c1e", "#17615e"))+
  xlab("")+
  ylab("Peak Variance")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

p_S

var_all_G2M <- var_all_long %>% filter(phase == "G2M")

stat_t.test_G2M <- var_all_G2M %>%
  t_test(value ~ FM_names, ref.group = "All genes in G2M phase")

p_G2M <- ggplot(var_all_G2M)+
  geom_errorbar(aes(x = reorder(FM_names, med, decreasing = T),
                    color = FM_names, ymin = low, ymax = high), width = 0.2)+
  geom_bar(aes(x = reorder(FM_names, med, decreasing = T),
               fill = FM_names, color = FM_names,
               y = med), data = unique(var_all_G2M[,c(1,5:7)]),
           width = 0.8, alpha = 0.5,
           stat = "identity")+
  geom_text(aes(reorder(FM_names, med, decreasing = T),
                y = 0.0002,
                label = format(med, scientific=TRUE, digits=2)),
            size = 2)+
  scale_y_continuous(position = "right")+
  geom_text(data = stat_t.test_G2M[match(levels(reorder(var_all_G2M$FM_names, var_all_G2M$med,
                                                      decreasing = T))[-1],
                                       stat_t.test_G2M$group2),],
            aes(x = 2:9, y = 0.010, label = p.adj.signif))+
  # 颜色模式：
  scale_fill_manual(values = c("#dd67a7", "#f0ad9c", "#85c2a3",
                               "#ffb4ec", "#cb773d",
                               "#aec221", "#70be6f", "#ef9c1e", "#17615e"))+
  scale_color_manual(values = c("#dd67a7", "#f0ad9c", "#85c2a3",
                                "#ffb4ec", "#cb773d",
                                "#aec221", "#70be6f", "#ef9c1e", "#17615e"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        # panel.background = element_rect(fill = "#fff4da"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

p_G2M


library(cowplot)

p_EF <- plot_grid(p_S, p_G2M, rel_widths = c(7, 9))
