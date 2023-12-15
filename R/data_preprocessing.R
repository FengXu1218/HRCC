## This script is a collection of functions used in data processing ---------------

#' remove_pc_outliers
#'
#' @param pc_score A dataframe containing scores on PC1 and PC2
#'
#' @return
#' @export
#'
#' @examples
remove_pc_outliers <- function(pc_score, time.iqr = 1){
  find_outlier_index <- function(x) {
    outlier.low <- quantile(x,probs=c(0.25)) - IQR(x)*time.iqr
    outlier.high <- quantile(x,probs=c(0.75)) + IQR(x)*time.iqr
    index <- which(x > outlier.high|x < outlier.low)
    return(index)
  }
  index1 <- find_outlier_index(pc_score[,1])
  index2 <- find_outlier_index(pc_score[,2])

  final_index <- setdiff(1:nrow(pc_score),union(index1, index2))
  return(final_index)
}

#' pca_plot
#'
#' This function is used to draw the distribution map of each point
#' on pc1 and pc2 after the first step of pca, and is used to screen
#' and remove outliers
#'
#' @param data Expression matrix of cell cycle genes
#' @param time.iqr times of IQR, defualt is 1.5
#'
#' @return
#' @export
#'
#' @examples
pca_plot <- function(data, time.iqr = 1){
  pca <- prcomp(t(data))
  score <- pca$x
  score <- as.data.frame(score[,1:2])
  final_index <- remove_pc_outliers(score, time.iqr = time.iqr)
  score$outliers <- ifelse(1:nrow(score) %in% final_index, "FALSE", "TRUE")

  # scatter plot of pca:
  p1 <- ggplot(score)+
    geom_point(aes(PC1, PC2, color = outliers))+
    scale_color_manual(values = c("TRUE" = "#fc011a", "FALSE" = "#75aadb"))+
    theme_bw()
  p1

  # Plot the distribution of points on pc1 and pc2:
  p2 <- ggplot(score)+
    geom_point(aes(PC1, 1:nrow(score), color = outliers))+
    scale_color_manual(values = c("TRUE" = "#fc011a", "FALSE" = "#75aadb"))+
    xlab("PC1")+
    ylab("index")+
    theme_bw()

  p3 <- ggplot(score)+
    geom_point(aes(PC2, 1:nrow(score), color = outliers))+
    scale_color_manual(values = c("TRUE" = "#fc011a", "FALSE" = "#75aadb"))+
    xlab("PC2")+
    ylab("index")+
    theme_bw()

  p <- p2/p3

  return(list(plot1 = p1, plot2 = p, score = score, final_index = final_index))
}


#' rank_to_angle
#'
#' @param rank_data a vector of cell rank
#' @param max_rank the length of cell
#'
#' @return
#' @export
#'
#' @examples
rank_to_angle <- function(rank_data, max_rank){
  angle_data <- rank_data*((2*pi)/max_rank) - pi
  return(angle_data)
}


#' counts2TPM
#'
#' @param count
#' @param efflength
#'
#' @return
#' @export
#'
#' @examples
counts2TPM <- function(count, efflength){
  RPK <- count/(efflength/1000)
  PMSC_rpk <- sum(RPK)/1e6
  TPM <- RPK/PMSC_rpk
  return(TPM)
}
