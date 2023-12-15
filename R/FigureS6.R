###### 差异分析 --------------
library(Seurat)
library(CAFF)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggplotify)
library(ClusterGVis)
library(org.Hs.eg.db)
library(clusterProfiler)

dir_path <- "DataSets/nolabel/GSE129447_HeLa1/"
load("Rdata/HeLa1_phase_correct.Rdata")

make_Trans_waves <- function(phase = 3, data = "HeLa1"){
  # 载入数据：
  load(paste0(dir_path, "Rdata/data_new.Rdata"))

  # 根据分期结果重新排序：
  data_cc_ordered <- data_cc_ordered[,names(HeLa1_phase_correct[[phase-2]])]

  # 重做傅里叶滤波：
  data_cc_bind <- cbind(data_cc_ordered, data_cc_ordered, data_cc_ordered)
  X_fftf <- data_cc_bind

  for (i in 1:nrow(data_cc_bind)) {
    res <- fftf(1:ncol(data_cc_bind), exp(as.numeric(data_cc_bind[i,])), 1/(ncol(data_cc_bind)/6), plot = F)
    X_fftf[i,] <- res[[1]]
    print(i)
  }

  HeLa1_fftf <- X_fftf[, ((ncol(X_fftf)/3+1):(ncol(X_fftf)*2/3))]

  ########## 周期蛋白的波动图 -----------
  tmp_data <- as.data.frame(HeLa1_fftf[c("CCNE1", "CCNB1", "CCNA2","PCNA"),])
  tmp_data$gene <- rownames(tmp_data)

  tmp_data <- tmp_data %>% pivot_longer(-gene,
                                        names_to = "cell",
                                        values_to = "value")
  tmp_data <- as.data.frame(tmp_data)

  tmp_data$cell <- factor(tmp_data$cell, levels = colnames(HeLa1_fftf))

  line_data <- cumsum(table(HeLa1_phase_correct[[phase-2]]))

  p <- ggplot(tmp_data)+
    geom_line(aes(cell, exp(value), color = gene, group = gene),
              linewidth = 0.6)+
    geom_vline(xintercept = line_data, linetype = "dashed")+
    scale_color_manual(values = c("#6883af", "#ef9973", "#cd6c72", "#6db97d"))+
    xlab("")+
    ylab("")+
    ggtitle(paste0("Phase = ", phase))+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = unit(6, "pt")),
          axis.ticks.x = element_blank(),
          axis.title = element_text(face = "bold.italic"),
          plot.title = element_text(hjust = 0.5, size = unit(8, "pt"),
                                    face = "bold.italic"),
          legend.position = "none")

  return(p)
}

p_list <- list()
for (phase in 3:10) {
  p_list[[phase-2]] <- make_Trans_waves(phase = phase)
}

library(cowplot)
p_all <- plot_grid(plotlist = p_list, ncol = 2)
p_all
