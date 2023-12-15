######### 三组基因的ridge plot ---------------
# 绘图：
load("Rdata/peak_genes_3Pattern.Rdata")
load("Rdata/max_index_norm_152.Rdata")

Peak_names <- ls(pattern = "Peaks")

for (j in 1:length(Peak_names)) {
  tmp <- max_index_norm[get(Peak_names[j]), ]
  tmp <- tmp[levels(reorder(get(Peak_names[j]), density_data[which(density_data$gene %in% get(Peak_names[j])),]$x)), ]
  p_list <- list()

  p <- ggplot()+
    geom_density(aes(na.omit(tmp[1,])), color = "#f064af",
                 fill = "#f064af", alpha= 0.6, bw = 0.03)+
    ylab(rownames(tmp)[1])+
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
    colnames(tmp2)[i] <- "Time series (t)"
    p_list[[i]] <- ggplot(tmp2)+
      geom_density(aes(x = `Pseudo Time`),
                   color = color[i],bw = 0.03,
                   fill = color[i], alpha= 0.6)+
      ylab(rownames(tmp)[i])+
      xlab("")+
      theme_minimal()+
      theme(panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5))+
      xlim(0,1)
  }

  p_all <- plot_grid(plotlist = p_list, ncol = 3)

  assign(paste0("p", Peak_names[j]), p_all)

  ggsave(plot = p_all, filename = paste0("Figures/supplymentary figures/p_", Peak_names[j], ".pdf"),
         height = nrow(tmp)*2/9, width = 10)
}


load("Rdata/peak_genes_3Pattern.Rdata")
load("Rdata/max_index_norm_152.Rdata")

Peak_names <- ls(pattern = "Peaks")

for (j in 1:length(Peak_names)) {
  tmp <- max_index_norm[get(Peak_names[j]), ]
  tmp <- tmp[levels(reorder(get(Peak_names[j]), density_data[which(density_data$gene %in% get(Peak_names[j])),]$x)), ]
  p_list <- list()

  p <- ggplot()+
    geom_density(aes(na.omit(tmp[1,])), color = "#f064af",
                 fill = "#f064af", alpha= 0.6, bw = 0.03)+
    # ylab(rownames(tmp)[1])+
    ylab("")+
    xlab("")+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5))+
    xlim(0,1)
  # p

  p_list[[1]] <- p

  color <- colorRampPalette(c("#f064af" , "#e4daa3", "#74c464"
                              ))(nrow(tmp))

  for (i in 2:nrow(tmp)) {
    tmp2 <- as.data.frame(t(tmp))
    colnames(tmp2)[i] <- "Time series (t)"
    p_list[[i]] <- ggplot(tmp2)+
      geom_density(aes(x = `Pseudo Time`),
                   color = color[i],bw = 0.03,
                   fill = color[i], alpha= 0.6)+
      # ylab(rownames(tmp)[i])+
      ylab("")+
      xlab("")+
      theme_minimal()+
      theme(panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5))+
      xlim(0,1)
  }

  p_all <- plot_grid(plotlist = p_list, ncol = 3)

  assign(paste0("p", Peak_names[j]), p_all)

  ggsave(plot = p_all, filename = paste0("Figures/supplymentary figures/p_nolabel_", Peak_names[j], ".pdf"),
         height = nrow(tmp)*2/9, width = 10)
}


