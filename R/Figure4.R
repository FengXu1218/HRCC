######################### Figure4 -----------------
load("Rdata/test_acc.Rdata")

########## p_A -------------
data_A <- data.frame(Model = c("XGBoost", "SVM", "RF", "NN", "Ensemble"),
                     Accuracy = c(rowMeans(base_model_acc_list[[7]]),
                                  test_confusion_mat[[7]][["overall"]][["Accuracy"]]))

p_A <- ggplot(data_A)+
  geom_linerange(aes(x = reorder(Model, Accuracy, decreasing = T),
                     ymin = 0, ymax = Accuracy), linetype = "dashed")+
  geom_point(aes(reorder(Model, Accuracy, decreasing = T),
                 Accuracy, fill = Model),
             size = 6, shape = 21)+
  geom_text(aes(reorder(Model, Accuracy, decreasing = T),
                Accuracy+0.1,
                label = paste0(round(Accuracy, 4)*100, "%")),
            size = 2)+
  scale_fill_brewer(palette = "Set2")+
  ylim(0, 1)+
  xlab("")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, face = "bold.italic"),
        axis.title.y = element_text(size = 8, face = "bold.italic"))
p_A

########### p_B&C -------------
# phase = 6:
pred_mat <- as.data.frame(test_confusion_mat[[4]][["table"]])
colnames(pred_mat) <- c("Predictions", "TICC", "Freq")

p_B <- ggplot(pred_mat, aes(TICC, Predictions))+
  geom_point(aes(size = Freq), color = "#83c6e8")+
  geom_text(aes(label = Freq), size = 3)+
  geom_hline(yintercept = 0.5:6.5, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0.5:6.5, linetype = "dashed", color = "grey")+
  scale_size_continuous(range = c(0, 8))+
  scale_x_discrete(expand = c(0.085, 0.085), labels = paste0("C6.", 1:6))+
  scale_y_discrete(expand = c(0.085, 0.085), labels = paste0("C6.", 1:6))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = "bold.italic", size = 8),
        axis.text = element_text(size = 6),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")
p_B

# phase = 9
pred_mat <- as.data.frame(test_confusion_mat[[7]][["table"]])
colnames(pred_mat) <- c("Predictions", "TICC", "Freq")

p_C <- ggplot(pred_mat, aes(TICC, Predictions))+
  geom_point(aes(size = Freq), color = "#83c6e8")+
  geom_text(aes(label = Freq), size = 3)+
  geom_hline(yintercept = 0.5:9.5, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0.5:9.5, linetype = "dashed", color = "grey")+
  scale_size_continuous(range = c(0, 8))+
  scale_x_discrete(expand = c(0.057, 0.057), labels = paste0("C9.", 1:9))+
  scale_y_discrete(expand = c(0.057, 0.057), labels = paste0("C9.", 1:9))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = "bold.italic", size = 8),
        axis.text = element_text(size = 6),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")
p_C



################ p_D & E ----------------------
library(tidyverse)
library(Seurat)
library(Matrix)
library(caret)
library(keras)
library(tensorflow)
library(pROC)
library(e1071)
library(class)
library(randomForest)
library(xgboost)
library(plotrix)

dir_path <- "DataSets/cell_cycle_arrest/GSE211617/"
load("Rdata/ensemble_model_9.Rdata")
load("Rdata/gene_list.Rdata")
load(paste0(dir_path, "Rdata/exp_cpm_cc.Rdata"))

source("R/get_meta_test_mat.R")

use_python("D:\\Program Files\\Anaconda3", required = T)

meta_test_matrix <- get_meta_test_mat(base_model_all, nn_model_dir = "Rdata/nn_models/Phase_9",
                                      data_test_meta = t(scale(all_test)))

Pred <- predict(meta_model, newdata = meta_test_matrix)

group_names <- rep(c("IR2h", "IR6h", "Ctrl"),
                   c(ncol(IR2h_cpm_cc)-1, ncol(IR6h_cpm_cc)-1, ncol(Ctrl_cpm_cc)-1))

Pred_IR2h <- Pred[group_names == "IR2h"]
Pred_IR6h <- Pred[group_names == "IR6h"]
Pred_Ctrl <- Pred[group_names == "Ctrl"]

tab_IR2h <- table(Pred_IR2h)
tab_IR2h <- tab_IR2h/sum(tab_IR2h)
tab_IR6h <- table(Pred_IR6h)
tab_IR6h <- tab_IR6h/sum(tab_IR6h)
tab_Ctrl <- table(Pred_Ctrl)
tab_Ctrl <- tab_Ctrl/sum(tab_Ctrl)

plot_data <- data.frame(Phase = rep(1:9, 3),
                        Propotion = c(tab_IR2h, tab_IR6h, tab_Ctrl),
                        group = rep(c("IR_2h", "IR_6h", "Ctrl"), each = 9))

plot_data$group <- factor(plot_data$group, levels = c("Ctrl", "IR_2h", "IR_6h"))

compute_std <- function(prop, n){
  std_p <- c()
  for (i in 1:length(prop)) {
    std_p[i] <- sqrt(prop[i]*(1-prop[i])/n)
  }

  return(std_p)
}

std_IR2h <- compute_std(prop = tab_IR2h, n = length(Pred_IR2h))
std_IR6h <- compute_std(prop = tab_IR6h, n = length(Pred_IR6h))
std_Ctrl <- compute_std(prop = tab_Ctrl, n = length(Pred_Ctrl))

plot_data$std <- c(std_IR2h, std_IR6h, std_Ctrl)

# barplot
load("Rdata/Ensemble_all_pred_IR.Rdata")

p_GSE211617 <- ggplot(plot_data)+
  geom_errorbar(aes(x = Phase, ymin = Propotion - std,
                    ymax = Propotion + std, color = group),
                position = position_dodge(0.9), width = 0.5)+
  geom_col(aes(Phase, Propotion, fill = group, color = group),
           alpha = 0.6, position = position_dodge(0.9), width = 0.8)+
  geom_text(aes(x = Phase, y = Propotion+0.02,
                label = paste0(round(Propotion, 4)*100, "%"), group = group),
            position = position_dodge(0.9), size = 2)+
  scale_color_manual(name = "Group", values = c("#da6f71", "#efd682","#75c078"))+
  scale_fill_manual(name = "Group", values = c("#da6f71", "#efd682","#75c078"))+
  scale_x_continuous(breaks = 1:9, labels = paste0("C9.", 1:9))+
  scale_y_continuous(breaks = seq(0, 0.3, 0.1),
                     labels = paste0(seq(0, 30, 10), "%"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.key.size = unit(10, "pt"),
        strip.text = element_text(face = "bold.italic"),
        title = element_text(face = "bold.italic", size = 8),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(face = "bold.italic"))
p_GSE211617


# line plot:
# Ctrl：
Ctrl_cpm_cc <- Ctrl_cpm_cc %>% column_to_rownames("geneName")
Ctrl_sce <- CreateSeuratObject(Ctrl_cpm_cc)
Ctrl_sce$cell.cycle.phase <- paste0("C9.", as.numeric(Pred_Ctrl))
Ctrl_sce <- ScaleData(Ctrl_sce)

Ctrl_feature <- as.data.frame(Ctrl_sce@assays[["RNA"]]@scale.data[c("CCNE1", "PCNA", "CCNA2", "CCNB1"), ]) %>%
  rownames_to_column("Gene_name") %>%
  pivot_longer(cols = !Gene_name, names_to = "cell_name", values_to = "expression")
Ctrl_feature$phase <- rep(as.character(Ctrl_sce$cell.cycle.phase), 4)

Ctrl_feature <- Ctrl_feature %>%
  group_by(Gene_name, phase) %>%
  mutate(mean_exp = mean(expression),
         group = "Control",
         sd = std.error(expression),
         phase = factor(phase, levels = paste0("C9.", 1:9)),
         Gene_name = factor(Gene_name, levels = c("CCNE1", "PCNA", "CCNA2", "CCNB1"))) %>%
  ungroup()

# IR2h:
IR2h_cpm_cc <- IR2h_cpm_cc %>% column_to_rownames("geneName")
IR2h_sce <- CreateSeuratObject(IR2h_cpm_cc)
IR2h_sce$cell.cycle.phase <- paste0("C9.", as.numeric(Pred_IR2h))
IR2h_sce <- ScaleData(IR2h_sce)

IR2h_feature <- as.data.frame(IR2h_sce@assays[["RNA"]]@scale.data[c("CCNE1", "PCNA", "CCNA2", "CCNB1"), ]) %>%
  rownames_to_column("Gene_name") %>%
  pivot_longer(cols = !Gene_name, names_to = "cell_name", values_to = "expression")
IR2h_feature$phase <- rep(as.character(IR2h_sce$cell.cycle.phase), 4)

IR2h_feature <- IR2h_feature %>%
  group_by(Gene_name, phase) %>%
  mutate(mean_exp = mean(expression),
         group = "IR-2h",
         sd = std.error(expression),
         phase = factor(phase, levels = paste0("C9.", 1:9)),
         Gene_name = factor(Gene_name, levels = c("CCNE1", "PCNA", "CCNA2", "CCNB1"))) %>%
  ungroup()


# IR6h 数据处理：
IR6h_cpm_cc <- IR6h_cpm_cc %>% column_to_rownames("geneName")
IR6h_sce <- CreateSeuratObject(IR6h_cpm_cc)
IR6h_sce$cell.cycle.phase <- paste0("C9.", as.numeric(Pred_IR6h))
IR6h_sce <- ScaleData(IR6h_sce)

IR6h_feature <- as.data.frame(IR6h_sce@assays[["RNA"]]@scale.data[c("CCNE1", "PCNA", "CCNA2", "CCNB1"), ]) %>%
  rownames_to_column("Gene_name") %>%
  pivot_longer(cols = !Gene_name, names_to = "cell_name", values_to = "expression")
IR6h_feature$phase <- rep(as.character(IR6h_sce$cell.cycle.phase), 4)

IR6h_feature <- IR6h_feature %>%
  group_by(Gene_name, phase) %>%
  mutate(mean_exp = mean(expression),
         group = "IR-6h",
         sd = std.error(expression),
         phase = factor(phase, levels = paste0("C9.", 1:9)),
         Gene_name = factor(Gene_name, levels = c("CCNE1", "PCNA", "CCNA2", "CCNB1"))) %>%
  ungroup()

plot_data <- rbind(Ctrl_feature, IR2h_feature, IR6h_feature)
plot_data$group <- factor(plot_data$group, levels = c("Control", "IR-2h", "IR-6h"))

p2 <- ggplot(unique(plot_data[,c("phase", "mean_exp",
                                 "Gene_name", "group")]),
             aes(phase, mean_exp)) +
  geom_line(aes(group = group, color = group), alpha = 0.8) +
  geom_point(aes(group = group, color = group))+
  scale_color_manual(values = c("#da6f71", "#efd682","#75c078"))+
  facet_wrap(Gene_name~., ncol = 4, strip.position = "top",
             scales = "free_y")+
  ylab(latex2exp::TeX("Normalized $log_{2}$(CPM)", bold = T, italic = T))+
  xlab("")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_blank(), #去除标题框背景
        strip.text = element_text(face = "bold.italic"),
        title = element_text(face = "bold.italic", size = 8),
        axis.text = element_text(size = 6),
        axis.text.x = element_text(face = "bold.italic",
                                   angle = 90, vjust = 0.5))

p2

library(patchwork)

p_GSE211617/p2 + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(2.5, 1))

