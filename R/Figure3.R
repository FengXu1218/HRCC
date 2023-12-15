###################### Figure 3 ####################
###### Fig 3A: heatmap of navid NMI ---------------
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)

source("R/visualization.R")

load("Rdata/phi_ticc.Rdata")

plot_data <- t(scale(temp_result[,!grepl("BC_SUM", colnames(temp_result))]))
# plot(3:20, apply(temp_result, 1, median))

plot_data <- plot_data[order(apply(plot_data, 1, which.max)), ]
colnames(plot_data) <- 3:20

col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("#74c464", "#fdedf6","#f064af"))

Heatmap(plot_data,
        # 设置颜色：
        col = col_fun,
        cluster_rows = F,
        cluster_columns = F,

        row_names_gp = gpar(fontsize = 5,
                            fontface = "italic" # 调整字体为斜体
        ),
        column_names_gp = gpar(fontsize = 5,
                               fontface = "italic" # 调整字体为斜体
        ),
        column_names_rot = 0,
        rect_gp = gpar(col = "white", lwd = 1),
        # 调整图例：
        show_heatmap_legend = FALSE,
        column_title = "Number of Stage",
        column_title_gp = gpar(fontsize = 8,
                               fontface = "bold.italic" # 调整字体为斜体
        ),
        column_title_side = "bottom")

###### Fig 3B&C: 6、9-phases TICC and FUCCI ---------------
## hESC:
load("Rdata/fftf_6.Rdata")
load("Rdata/stage_all_correct.Rdata")

fucci_hESC <- str_split(colnames(hESC_fftf_6), "_", simplify = T)[,1]
ticc_hESC <- sort(hESC_stage$cluster6)

mat_data_hESC <- data.frame(FUCCI = fucci_hESC, TICC = ticc_hESC)
mat_data_hESC <- mat_data_hESC %>% filter(FUCCI != "H1")

plot_data_hESC <- as.data.frame(table(mat_data_hESC))
plot_data_hESC$FUCCI <- factor(plot_data_hESC$FUCCI, levels = c("G2", "S", "G1"))
plot_data_hESC$TICC <- factor(plot_data_hESC$TICC, levels = c("F", "A", "B", "C", "D", "E"))

mat_point_plots(plot_data_hESC, x = "TICC", y = "FUCCI", phase = 6)


load("Rdata/fftf_9.Rdata")

fucci_hESC <- str_split(colnames(hESC_fftf_9), "_", simplify = T)[,1]
ticc_hESC <- sort(hESC_stage$cluster9)

mat_data_hESC <- data.frame(FUCCI = fucci_hESC, TICC = ticc_hESC)
mat_data_hESC <- mat_data_hESC %>% filter(FUCCI != "H1")

plot_data_hESC <- as.data.frame(table(mat_data_hESC))
plot_data_hESC$FUCCI <- factor(plot_data_hESC$FUCCI, levels = c("G2", "S", "G1"))
plot_data_hESC$TICC <- factor(plot_data_hESC$TICC, levels = c("G", "H", "I",
                                                              LETTERS[1:6]))

mat_point_plots(plot_data_hESC, x = "TICC", y = "FUCCI", phase = 9,
                expand_x = c(0.058, 0.06))
