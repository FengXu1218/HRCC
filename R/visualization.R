#' pca_plot2
#'
#'
#' This function is used to plot a scatter plot of the expression matrix distribution on PC1 and PC2
#' @param data expression matrix
#'
#' @return
#' @export
#'
#' @examples
pca_plot2 <- function(data){
  pca <- prcomp(t(data))
  score <- pca$x

  p <- ggplot(as.data.frame(score))+
    geom_point(aes(PC1, PC2), color = "red")+
    theme_bw()

  x <- score[,1:2]

  return(list(plot = p, x = x))
}


#' pca_scatter_plot
#'
#'
#' make a pca scatter plot
#'
#' @param pca_list the list after pca by function prcomp.
#' @param title title
#' @param sub_title sub title
#'
#' @return
#' @export
#'
#' @examples
pca_scatter_plot <- function(pca_list, title, sub_title = "",
                             add_legend = T,
                             colors = c("#3d4ba0", "#bb6bab", "#e1cf5b", "#499a44","#3d4ba0")){
  data <- as.data.frame(pca_list$x)

  # convert rectangular coordinates to polar coordinates:
  THETA <- apply(data[,1:2], 1, cart2pol)[1,]
  data$THETA <- THETA

  pc_var <- summary(pca_list)$importance
  pc1_var <- paste0("PC1(",round(pc_var[2,1]*100,2),"%)")
  pc2_var <- paste0("PC2(",round(pc_var[2,2]*100,2),"%)")
  p1 <- ggplot(data, aes(PC1, PC2))+
    geom_point(aes(color = THETA))+
    scale_color_gradientn(colors = colors,
                          # limits = c(-pi, pi),
                          # breaks = c(-pi, -0.5*pi, 0, 0.5*pi, pi),
                          # labels = expression(0, pi/2, pi,
                          #                     3*pi/2, 2*pi)
    ) +
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold.italic"),
      plot.subtitle = element_text(hjust = 0.5, face = "bold.italic"),
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold.italic"),
      legend.position = "none"
    )+
    xlab(pc1_var)+
    ylab(pc2_var)+
    ggtitle(title, subtitle = sub_title)

  if (add_legend == TRUE) {
    legend_data <- data.frame(x = c(rep(0, 100), rep(1,100)), y = c(rep(0, 100), rep(1,100)),
                              group = rep(1:100, 2))
    p2 <- ggplot(legend_data)+
      geom_col(aes(x, y, fill = group))+
      scale_fill_gradientn(colors = colors)+
      scale_y_continuous(breaks = seq(0,100,25),
                         labels = expression(0, pi/2, pi,
                                             3*pi/2, 2*pi))+
      coord_polar(theta = "y", start = pi/2, direction = -1)+
      theme_minimal()+
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 5),
            axis.text.y = element_blank(),
            legend.position = "none")+
      xlab("")+
      ylab("")

    library(cowplot)

    p <- ggdraw() +
      draw_plot(p1, 0, 0, 0.8, 1) +
      draw_plot(p2, 0.75, 0.6, 0.3, 0.3)
    return(p)
  }

  return(p1)
}


#' scatter_plot
#'
#' This function is used to make the scatter plot of the expression
#' varying with the pseudotime.
#'
#' @param exp_mat expression matrix
#' @param gene_list a list of genes which you want to display
#' @param n_col column number of the plot
#' @param y_lab the label of y-axis
#'
#' @return
#' @export
#'
#' @examples
scatter_plot <- function (exp_mat, gene_list, y_lab = "Variance", n_col = 5)
{
  library(tidyverse)
  if (length(gene_list) == 1) {
    data_visual <- as.data.frame(t(exp_mat[gene_list, ]))
    rownames(data_visual) <- gene_list
  } else {
    data_visual <- as.data.frame(exp_mat[gene_list, ])
  }
  data_visual$gene <- rownames(data_visual)
  data_visual_long <- pivot_longer(data_visual, cols = !gene,
                                   names_to = "Order", values_to = "Variance")
  data_visual_long$Order <- rep(1:(ncol(data_visual) - 1),
                                nrow(data_visual))
  p_list <- list()
  for (i in 1:nrow(data_visual)) {
    gene <- rownames(data_visual)[i]
    p <- ggplot(data_visual_long[data_visual_long$gene ==
                                   gene, ], aes(Order, Variance)) + geom_point(color = "#988d7b") +
      geom_smooth(method = "loess", formula = "y ~ x",
                  fill = "#b2e7fa", color = "#00aeef", alpha = 0.6,
                  se = 0.98) + ylab(y_lab) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                panel.grid = element_blank(), axis.title = element_text(face = "bold.italic"),
                  ) + ggtitle(gene)
    p_list[[i]] <- p
  }
  library(cowplot)
  p <- plot_grid(plotlist = p_list, ncol = n_col)
  return(list(data_visual = data_visual, data_visual_long = data_visual_long,
              plot = p))
}





#' heatmap_plot
#'
#'
#' This function is used to reorder the Fourier filtered results
#' and make a heatmap.
#'
#' @param X_fftf2 Matrix after fftf
#' @param name title of the heatmap
#' @param output weather to output
#' @param save_path if output = TRUE, you need to provide the save path
#'
#' @return
#' @export
#'
#' @examples
heatmap_plot <- function(X_fftf2, name,
                         output = FALSE,
                         save_path = "."){
  final_index <- get_row_order(X_fftf2)[[2]]

  p <- pheatmap::pheatmap(as.matrix(X_fftf2[final_index,]),
                          # color = colorRampPalette(c("#14315c","#9ec4dc", "#ffffff", "#c66650", "#5f101f"))(100),
                          cluster_rows = F,
                          fontsize_row = 5,
                          show_colnames = F,
                          cluster_cols = F,
                          main = paste("Reconstructed signal of", name),
                          legend = F)

  if (output == T) {
    pdf(paste0(save_path, name, "_heatmap.pdf"), height = 7, width = 5)
    print(p)
    dev.off()
  }
  return(list(plot = p, final_index = final_index))
}


#' cor_scatter_plot
#'
#'
#' This function is used to plot a scatter plot that maximizes correlation
#'
#' @param index_data a data frame.
#'
#' @return
#' @export
#'
#' @examples
cor_scatter_plot <- function(index_data){
  p_list <- list()
  library(Hmisc)

  t <- 1
  for (i in 1:(ncol(index_data)-1)) {
    for (j in 2:ncol(index_data)) {
      if (j > i) {
        res_out <- rings2aligment(index_data[,j], index_data[,i])
        res <- rcorr(index_data[,j][res_out$output_posi], index_data[,i])
        p_value <- signif(res$P[1,2], 2)
        cor_value <- round(res$r[1,2], 2)

        index_data_new <- index_data[,c(j,i)]
        index_data_new[,1] <- index_data_new[res_out$output_posi,1]
        colnames(index_data_new) <- c("x", "y")

        p <- ggplot(index_data_new, aes(x = x, y = y))+
          geom_point(color = "#988d7b")+
          geom_smooth(method = "lm", formula = y ~ x,
                      fill = "#b2e7fa", color = "#00aeef", alpha = 0.8)+
          theme_bw()+
          xlab(colnames(index_data)[j])+
          ylab(colnames(index_data)[i])+
          theme(
            panel.grid = element_blank(),
            axis.title = element_text(face = "bold.italic"),
            plot.title = element_text(hjust = 0.5, size = 10)
          )+
          labs(title = paste0("r =", cor_value, ", q = ", p_value))
        p_list[[t]] <- p
        t <- t+1
      }
    }
  }

  return(p_list)
}



#' circ_cor_scatter_plot
#'
#'
#' This function is used to draw a scatter plot of circular correlation coefficient for circle data
#'
#' @param circ_data a data frame with colnames
#' @param cor_method method of calculating correlation
#'
#' @return
#' @export
#'
#' @examples
circ_cor_scatter_plot <- function (circ_data, cor_method = c("circ_cor", "Spearman_cor"))
{
  p_list <- list()
  data_list <- list()
  cor_mat <- matrix(NA, nrow = ncol(circ_data), ncol = ncol(circ_data))
  rownames(cor_mat) <- colnames(circ_data)
  colnames(cor_mat) <- colnames(circ_data)
  p_mat <- matrix(NA, nrow = ncol(circ_data), ncol = ncol(circ_data))
  rownames(p_mat) <- colnames(circ_data)
  colnames(p_mat) <- colnames(circ_data)
  library(BAMBI)
  t <- 1
  for (i in 1:(ncol(circ_data) - 1)) {
    for (j in (i+1):ncol(circ_data)) {
      circ_data_new <- circ_data[, c(i, j)]
      circ_data_new <- circ_data_new[order(circ_data_new[, 1], decreasing = T), ]
      for (n in 1:ncol(circ_data_new)) {
        circ_data_new[, n] <- rank_to_angle(circ_data_new[, n],
                                            max_rank = max(circ_data_new[, n]))
      }
      res_out <- rings2aligment(circ_data_new[, 1],
                                circ_data_new[, 2])
      circ_data_new[, 1] <- circ_data_new[, 1][res_out$output_posi]
      theta1 <- circ_data_new[, 1]
      theta2 <- circ_data_new[, 2]
      for (n in 1:nrow(circ_data)) {
        if (abs(theta1[n] - theta2[n]) > pi) {
          if (theta1[n] < theta2[n]) {
            theta1[n] <- 2 * pi + theta1[n]
          }
          else {
            theta2[n] <- 2 * pi + theta2[n]
          }
        }
      }
      circ_data_new <- as.data.frame(cbind(theta1,
                                           theta2))
      res_out <- rings2aligment(circ_data_new[, 1],
                                circ_data_new[, 2])
      circ_data_new[, 1] <- circ_data_new[, 1][res_out$output_posi]
      colnames(circ_data_new) <- c("x", "y")
      if (cor_method == "circ_cor") {
        res <- circ_cor(as.matrix(circ_data_new),
                        type = "js")
        p_value <- round(attributes(res)$pval, 2)
        cor_value <- round(res, 2)
        cor_mat[j, i] <- cor_value
        cor_mat[i, j] <- cor_value
        p_mat[j, i] <- p_value
        p_mat[i, j] <- p_value
      }
      if (cor_method == "Spearman_cor") {
        res <- rcorr(circ_data_new[, 1], circ_data_new[,
                                                       2])
        p_value <- signif(res$P[1, 2], 2)
        cor_value <- round(res$r[1, 2], 2)
        cor_mat[j, i] <- cor_value
        cor_mat[i, j] <- cor_value
        p_mat[j, i] <- p_value
        p_mat[i, j] <- p_value
      }
      p <- ggplot(circ_data_new, aes(x = x, y = y)) +
        geom_point(color = "#988d7b") + geom_smooth(method = "lm",
                                                    formula = y ~ x, fill = "#b2e7fa", color = "#00aeef",
                                                    alpha = 0.8) + theme_bw() + xlab(colnames(circ_data)[j]) +
        ylab(colnames(circ_data)[i]) + theme(panel.grid = element_blank(),
                                             axis.title = element_text(face = "bold.italic"),
                                             plot.title = element_text(hjust = 0.5, size = 10)) +
        labs(title = paste0(cor_method, "=", cor_value,
                            ", p_value = ", p_value))
      p_list[[t]] <- p
      data_list[[t]] <- circ_data_new
      t <- t + 1
    }
  }
  diag(cor_mat) <- 1
  p_mat[p_mat == 0] <- 1e-15
  diag(p_mat) <- 0
  return(list(plot = p_list, data = data_list, cor_mat = cor_mat,
              p_mat = p_mat))
}




#' venn_upset_plot
#'
#' @param gene_list a list of genes.
#' @param save_path path to save the venn and upset plot
#' @param venn_name name of venn plot
#' @param upset_name name of upset plot
#'
#' @return
#' @export
#'
#' @examples
venn_upset_plot <- function(gene_list, save_path,
                            venn_name = "cc_genes_venn",
                            upset_name = "cc_genes_upset"){
  library(RColorBrewer)
  library(VennDiagram)

  category.names <- names(gene_list)
  fill_color <- brewer.pal(length(gene_list), "Set1")

  venn.diagram(
    x = gene_list,
    category.names = category.names,
    filename = paste0(save_path, venn_name, ".png"),

    imagetype="png" ,
    height = 1000 ,
    width = 1000 ,
    resolution = 300,
    compression = "lzw",

    col = "white",
    lty = 1,
    lwd = 1,
    fill = fill_color,
    alpha = 0.90,

    label.col = "black",
    cex = .5,
    fontfamily = "serif",
    fontface = "bold",

    cat.col = fill_color,
    cat.cex = .6,
    cat.fontfamily = "serif",

    disable.logging = T
  )
  message(paste0("Venn plot was save in ", save_path))

  # Upset-Plot
  library(UpSetR)

  # Plot
  pdf(paste0(save_path, upset_name, ".pdf"), height = 5, width = 10)
  print(upset(fromList(gene_list), order.by = "freq"))
  dev.off()

  message(paste0("upset plot was save in", save_path))
}


#' mat_point_plots
#'
#' @param data A data frame contained vars of x, y and Freq
#' @param x the column name of x
#' @param y the column name of y
#' @param size the column name of point size
#'
#' @return
#' @export
#'
#' @examples
mat_point_plots <- function(data, x, y, size = "Freq",
                            phase_number = 5,
                            expand_x = c(0.084, 0.09),
                            expand_y = c(0.17, 0.17)){
  colnames(data)[which(colnames(data) == x)] <- "x"
  colnames(data)[which(colnames(data) == y)] <- "y"

  p <- ggplot(data, aes(x, y))+
    geom_point(aes(size = Freq), color = "#81c3e2")+
    geom_text(aes(label = Freq), color = "black")+
    geom_vline(xintercept =seq(0.5, phase_number+0.5, 1), color = "#bbbbbb", linetype = "dashed")+
    geom_hline(yintercept=seq(0.5, 3.5, 1), color = "#bbbbbb", linetype = "dashed")+
    scale_size_continuous(range = c(0,10)) +
    scale_x_discrete(expand = expand_x)+
    scale_y_discrete(expand = expand_y)+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          legend.position = "none")
  p

  return(p)
}


#' mat_point_bar_plots
#'
#' @param data A data frame contained vars of x, y and Freq
#' @param x the column name of x
#' @param y the column name of y
#'
#' @return
#' @export
#'
#' @examples
mat_point_bar_plots <- function (data, x, y, expand_x = c(0.084, 0.09),
                                 expand_y = c(0.17, 0.17),
                                 phase_number = 5,
                                 x_lab = "",
                                 point_size = 0.4, point_stroke = 0.05,
                                 fills = c("#7391c9", "#6fbe44", "#eb3d67"))
{
  colnames(data)[which(colnames(data) == x)] <- "x"
  colnames(data)[which(colnames(data) == y)] <- "y"
  data_new <- data[, c("x", "y")]
  n = 1
  for (i in 1:nrow(data)) {
    if (data$Freq[i] == 0) {
      next
    }
    else {
      data_new[n:(n + data$Freq[i] - 1), ] <- data[i,
                                                   2:1]
      n = n + data$Freq[i]
    }
  }
  p1 <- ggplot(data_new, aes(x, y)) +
    geom_jitter(aes(fill = y), color = "white", size = point_size,
                stroke = point_stroke, alpha = 0.8, width = 0.4, height = 0.4,shape = 21) +
    geom_vline(xintercept = seq(0.5, phase_number+0.5, 1),
               size = unit(0.5, "pt"),
               color = "black", linetype = "dashed") +
    geom_hline(yintercept = seq(0.5, 3.5, 1),
               size = unit(0.5, "pt"),
               color = "black", linetype = "dashed") +
    scale_fill_manual(values = fills) +
    scale_x_discrete(expand = expand_x, labels = paste0("C", phase_number, ".", 1:phase_number)) +
    scale_y_discrete(expand = expand_y) +
    ylab("FUCCI")+
    xlab(paste0(phase_number, "-phases of TICC with ", x_lab))+
    theme_minimal() + theme(panel.grid = element_blank(),
                            axis.text = element_text(size = 5),
                            axis.title = element_text(size = 6, face = "bold.italic"),
                            legend.position = "none")

  p2 <- ggplot(data) + geom_col(aes(x, Freq, fill = y), position = "dodge") +
    geom_text(aes(x, Freq + 30, color = y, label = Freq),
              position = position_dodge(width = 0.9), size = 2) +
    scale_x_discrete(expand = c(0.1, 0.1)) + scale_y_continuous(limits = c(0,
                                                                           max(data$Freq) + 50), expand = c(0, 0.5)) + scale_fill_manual(name = "",
                                                                                                                                         values = fills) + scale_color_manual(name = "",
                                                                                                                                                                              values = c("#7391c9", "#6fbe44", "#eb3d67")) + ylab("Cell Number") +
    theme_classic() + theme(panel.grid = element_blank(),
                            axis.text = element_text(size = 5),
                            axis.title = element_text(size = 6, face = "bold.italic"),
                            axis.text.x = element_blank(),
                            axis.title.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            legend.position = "none")
  library(patchwork)
  p <- p2/p1 + plot_layout(heights = c(0.8, 2))
  return(p)
}


