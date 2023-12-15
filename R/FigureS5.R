############ Figure S5 ---------------
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)

source("R/visualization.R")
load("Rdata/fftf_6.Rdata")
load("Rdata/fftf_9.Rdata")
load("Rdata/stage_all_correct.Rdata")

######## Figure S5A: phase = 6 ----------------
fucci_hU2OS <- read.csv("DataSets/label/GSE146773_hU2OS/rawdata/phase.csv")
ticc_hU2OS <- hU2OS_stage[order(hU2OS_stage$cluster6),c(1,5)]

mat_data_hU2OS <- merge(ticc_hU2OS, fucci_hU2OS, by.x = "cell_name", by.y = "cell")[,-1][,2:1]
colnames(mat_data_hU2OS) <- c("FUCCI","TICC")
mat_data_hU2OS <- mat_data_hU2OS %>% filter(FUCCI != "NotAssigned")

plot_data_hU2OS <- as.data.frame(table(mat_data_hU2OS))
plot_data_hU2OS$FUCCI <- factor(plot_data_hU2OS$FUCCI, levels = c("G2M", "S-ph", "G1"))
plot_data_hU2OS$TICC <- factor(plot_data_hU2OS$TICC, levels = c("F", "A", "B", "C", "D", "E"))

mat_point_plots(plot_data_hU2OS, x = "TICC", y = "FUCCI", phase = 6)

######## Figure S5B: phase = 9 ----------------
ticc_hU2OS <- hU2OS_stage[order(hU2OS_stage$cluster9),c(1,8)]

mat_data_hU2OS <- merge(ticc_hU2OS, fucci_hU2OS, by.x = "cell_name", by.y = "cell")[,-1][,2:1]
colnames(mat_data_hU2OS) <- c("FUCCI","TICC")
mat_data_hU2OS <- mat_data_hU2OS %>% filter(FUCCI != "NotAssigned")

plot_data_hU2OS <- as.data.frame(table(mat_data_hU2OS))
plot_data_hU2OS$FUCCI <- factor(plot_data_hU2OS$FUCCI, levels = c("G2M", "S-ph", "G1"))
plot_data_hU2OS$TICC <- factor(plot_data_hU2OS$TICC, levels = c("G", "H", "I",
                                                                LETTERS[1:6]))
mat_point_plots(plot_data_hU2OS, x = "TICC", y = "FUCCI", phase = 9,
                expand_x = c(0.058, 0.06))

######## Figure S5C ----------------
# 3 - 6
tmp_data1 <- HeLa1_stage[,c(2,5)]

for (i in 1:7) {
  tmp_data1 <- apply(tmp_data1, 2, function(x) gsub(LETTERS[i], i, x))
}

tmp_data1 <- t(apply(tmp_data1, 1, as.numeric))
tmp_data1[,1] <- tmp_data1[,1] - 1
tmp_data1[,2] <- tmp_data1[,2] - 1 + 3

# 6 - 9
tmp_data2 <- HeLa1_stage[,c(5,8)]

for (i in 1:9) {
  tmp_data2 <- apply(tmp_data2, 2, function(x) gsub(LETTERS[i], i, x))
}
tmp_data2 <- t(apply(tmp_data2, 1, as.numeric))
tmp_data2[,1] <- tmp_data2[,1] - 1 + 3
tmp_data2[,2] <- tmp_data2[,2] - 1 + 9

plot_data <- as.data.frame(rbind(unique(tmp_data1), unique(tmp_data2)))
plot_data <- arrange(plot_data, V1, V2)

freq_data <- rbind(as.data.frame(table(apply(tmp_data1, 1,
                                             function(x) paste(x, collapse = "")))),
                   as.data.frame(table(apply(tmp_data2, 1,
                                             function(x) paste(x, collapse = "")))))

plot_data$Freq <- freq_data$Freq

links <- plot_data
colnames(links) <- c("source", "target", "value")

nodes <- data.frame(name = c(paste(3, 1:3, sep = "_"),
                             paste(6, 1:6, sep = "_"),
                             paste(9, 1:9, sep = "_")))

library(networkD3)

sankeyNetwork(Links = links, Nodes = nodes,
              # 指定source、target、value和name：
              Source = "source",
              Target = "target",
              Value = "value",
              NodeID = "name", # 节点的名字
              # 调整配置：
              fontSize = 12, # 节点的字体大小
              nodeWidth = 30, # 节点的宽度
              nodePadding = 8 # 节点之间的距离
)



