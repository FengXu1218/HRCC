################### Figuer S8 --------------
dir_path <- "DataSets/cell_cycle_arrest/GSE227555/"
load("Rdata/gene_list.Rdata")
load("Rdata/ensemble_model_9.Rdata")
load(paste0(dir_path, "Rdata/exp_cpm.Rdata"))

UT_cpm_cc <- data_UT_cpm[intersect(gene_8, rownames(data_UT_cpm)), ]
dim(UT_cpm_cc)

H50_cpm_cc <- data_50_cpm[intersect(gene_8, rownames(data_50_cpm)), ]
dim(H50_cpm_cc)

tmp_data <- as.data.frame(matrix(NA, nrow = ncol(data_test_meta), ncol = 1))
rownames(tmp_data) <- colnames(data_test_meta)
tmp_data <- tmp_data %>% rownames_to_column("geneName")
UT_cpm_cc <- as.data.frame(UT_cpm_cc) %>% rownames_to_column("geneName")
H50_cpm_cc <- as.data.frame(H50_cpm_cc) %>% rownames_to_column("geneName")

all_test <- Reduce(function(x, y) left_join(x, y, by = "geneName"),
                   list(tmp_data, UT_cpm_cc, H50_cpm_cc))
all_test <- all_test[,-2] %>% column_to_rownames("geneName")
table(is.na(all_test))
all_test[is.na(all_test)] <- 0
all_test <- as.matrix(all_test)

source("R/get_meta_test_mat.R")
meta_test_matrix <- get_meta_test_mat(base_model_all, nn_model_dir = "Rdata/nn_models/Phase_9",
                                      data_test_meta = t(scale(all_test)))

Pred <- predict(meta_model, newdata = meta_test_matrix)

group_names <- rep(c("Un-Treated", "H2O2-Treated"),
                   c(ncol(UT_cpm_cc)-1, ncol(H50_cpm_cc)-1))

Pred_UT <- Pred[group_names == "Un-Treated"]
Pred_H50 <- Pred[group_names == "H2O2-Treated"]

tab_UT <- table(Pred_UT)
tab_UT <- tab_UT/sum(tab_UT)
tab_H50 <- table(Pred_H50)
tab_H50 <- tab_H50/sum(tab_H50)

plot_data <- data.frame(Phase = rep(1:9, 2),
                        Propotion = c(tab_UT, tab_H50),
                        group = rep(c("Un-Treated", "H2O2-Treated"), each = 9))

plot_data$group <- factor(plot_data$group, levels = c("Un-Treated", "H2O2-Treated"))

compute_std <- function(prop, n){
  std_p <- c()
  for (i in 1:length(prop)) {
    std_p[i] <- sqrt(prop[i]*(1-prop[i])/n)
  }

  return(std_p)
}

std_UT <- compute_std(prop = tab_UT, n = length(Pred_UT))
std_H50 <- compute_std(prop = tab_H50, n = length(Pred_H50))

plot_data$std <- c(std_UT, std_H50)

p_GSE227555 <- ggplot(plot_data)+
  geom_errorbar(aes(x = Phase, ymin = Propotion - std,
                    ymax = Propotion + std, color = group),
                position = position_dodge(0.9), width = 0.5)+
  geom_col(aes(Phase, Propotion, fill = group, color = group),
           alpha = 0.6, position = position_dodge(0.9), width = 0.8)+
  geom_text(aes(x = Phase, y = Propotion+0.02,
                label = paste0(round(Propotion, 4)*100, "%"), group = group),
            position = position_dodge(0.9), size = 2)+
  scale_color_manual(name = "Group", values = c("#da6f71", "#75c078"))+
  scale_fill_manual(name = "Group", values = c("#da6f71", "#75c078"))+
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
p_GSE227555


# line plot
# UT:
UT_cpm_cc <- UT_cpm_cc %>% column_to_rownames("geneName")
UT_sce <- CreateSeuratObject(UT_cpm_cc)
UT_sce$cell.cycle.phase <- paste0("C9.", as.numeric(Pred_UT))
UT_sce <- ScaleData(UT_sce)

UT_feature <- as.data.frame(UT_sce@assays[["RNA"]]@scale.data[c("CCNE1", "PCNA", "CCNA2", "CCNB1"), ]) %>%
  rownames_to_column("Gene_name") %>%
  pivot_longer(cols = !Gene_name, names_to = "cell_name", values_to = "expression")
UT_feature$phase <- rep(as.character(UT_sce$cell.cycle.phase), 4)

UT_feature <- UT_feature %>%
  group_by(Gene_name, phase) %>%
  mutate(mean_exp = mean(expression),
         group = "Un-Treated",
         sd = std.error(expression),
         phase = factor(phase, levels = paste0("C9.", 1:9)),
         Gene_name = factor(Gene_name, levels = c("CCNE1", "PCNA", "CCNA2", "CCNB1"))) %>%
  ungroup()

# H50：
H50_cpm_cc <- H50_cpm_cc %>% column_to_rownames("geneName")
H50_sce <- CreateSeuratObject(H50_cpm_cc)
H50_sce$cell.cycle.phase <- paste0("C9.", as.numeric(Pred_H50))
H50_sce <- ScaleData(H50_sce)

H50_feature <- as.data.frame(H50_sce@assays[["RNA"]]@scale.data[c("CCNE1", "PCNA", "CCNA2", "CCNB1"), ]) %>%
  rownames_to_column("Gene_name") %>%
  pivot_longer(cols = !Gene_name, names_to = "cell_name", values_to = "expression")
H50_feature$phase <- rep(as.character(H50_sce$cell.cycle.phase), 4)

H50_feature <- H50_feature %>%
  group_by(Gene_name, phase) %>%
  mutate(mean_exp = mean(expression),
         group = "H2O2-Treated",
         sd = std.error(expression),
         phase = factor(phase, levels = paste0("C9.", 1:9)),
         Gene_name = factor(Gene_name, levels = c("CCNE1", "PCNA", "CCNA2", "CCNB1"))) %>%
  ungroup()

plot_data <- rbind(UT_feature, H50_feature)
plot_data$group <- factor(plot_data$group, levels = c("Un-Treated", "H2O2-Treated"))
p2 <- ggplot(unique(plot_data[,c("phase", "mean_exp",
                                 "Gene_name", "group")]),
             aes(phase, mean_exp)) +
  geom_line(aes(group = group, color = group), alpha = 0.8) +
  geom_point(aes(group = group, color = group))+
  scale_color_manual(values = c("#da6f71", "#75c078"))+
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
p_GSE227555/p2 + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(2.5, 1))


