###################### Figure 2全部图形 ####################
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplotify)
library(cowplot)

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

HeLa1_raw_scale <- t(scale(t(HeLa1_raw)))
max(HeLa1_raw_scale)
min(HeLa1_raw_scale)
col_fun <- colorRamp2(c(-1, -0.5, 0.5, 2), c("#74c464", "#badfb7", "#fdedf6","#f064af"))

p_B <- Heatmap(HeLa1_raw_scale,
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
p_B
