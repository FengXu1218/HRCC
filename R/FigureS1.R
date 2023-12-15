######## Figure S1 -----------
library(tidyverse)
library(CAFF)
library(cowplot)

load("Rdata/cc_genes.Rdata")
load("Rdata/gene_list.Rdata")
exp_files <- list.files("Rdata/All_Exp/",pattern = "Rdata")

for (i in 1:length(exp_files)) {
  exp_file <- exp_files[i]
  load(paste0("Rdata/All_Exp/", exp_file))
  if ("exp_cpm" %in% ls()) {
    assign(gsub("GSE[0-9]*_", "",
                str_split(exp_file, "[.]", simplify = T)[1]),
           exp_cpm)
    rm("exp_cpm", "exp")
  } else {
    assign(gsub("GSE[0-9]*_", "",
                str_split(exp_file, "[.]", simplify = T)[1]),
           exp_tpm)
    rm("exp_tpm", "exp")
  }
}

exp_vars <- ls(pattern = "_exp")[c(6:9,5,13,11,3:4,10,14,12)]
p_list <- list()
title_names <- c("HeLaCCL2_2019a", "HeLaCCL2_2019b", "HeLaS3_2020KO",
                 "HeLaS3_2020a", "HeLaS3_2020b", "hNPC_2021", "hFib_2022",
                 "hESC_2019a", "hESC_2019b", "hESC_2015",
                 "hU2OS_2020", "hiPSC_2019")

for (i in 1:12) {
  print(i)
  exp_tmp <- get(exp_vars[i])

  # pca scatter plot:
  p1 <- pca_scatter_plot(prcomp(t(exp_tmp[intersect(rownames(exp_tmp),
                                                    seurat_cc_genes),])),
                         add_legend = F,
                         title = title_names[i],
                         sub_title = "Seurat")

  p2 <- pca_scatter_plot(prcomp(t(exp_tmp[intersect(rownames(exp_tmp),
                                                    GO0007049_genes),])),
                         add_legend = F,
                         title = title_names[i],
                         sub_title = "GO0007049")

  p3 <- pca_scatter_plot(prcomp(t(exp_tmp[intersect(rownames(exp_tmp),
                                                    cyclebase3.0_genes),])),
                         add_legend = F,
                         title = title_names[i],
                         sub_title = "Cyclebase3.0")

  p4 <- pca_scatter_plot(prcomp(t(exp_tmp[intersect(rownames(exp_tmp),
                                                    gene_9),])),
                         add_legend = F,
                         title = title_names[i],
                         sub_title = "PFM")

  p_list[[i]] <- plot_grid(p1, p2, p3, p4, ncol = 4)
}

p <- plot_grid(plotlist = p_list, ncol = 2)

