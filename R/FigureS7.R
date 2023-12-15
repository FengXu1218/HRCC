############### Figure S7  -----------------------
load("Rdata/test_acc.Rdata")

# A.
data_A <- list()
for (i in 1:length(test_confusion_mat)) {
  data <- data.frame(Model = c("XGBoost", "SVM", "RF", "NN", "Ensemble"),
                     Accuracy = c(rowMeans(base_model_acc_list[[i]]),
                                  test_confusion_mat[[i]][["overall"]][["Accuracy"]]))
  # data$Accuracy[4] <-data$Accuracy[3]
  # data$Accuracy[3] <- data$Accuracy[2]

  data_A[[i]] <- data
}

data_A <- Reduce(rbind, data_A)
data_A$phase <- rep(3:10, each = 5)

data_A <- data_A %>%
  group_by(phase) %>%
  arrange(phase, desc(Accuracy)) %>%
  ungroup()

data_A$Model <- factor(data_A$Model,
                       levels = c("Ensemble", "XGBoost", "NN",
                                  "SVM", "RF"))
p_A <- ggplot(data_A)+
  annotate("rect", xmin = seq(3.5, 9.5, 2), xmax = seq(4.5, 10.5, 2),
           ymin = 0, ymax = 1, fill = "#ffd92f", alpha = 0.1)+
  annotate("rect", xmin = seq(2.5, 8.5, 2), xmax = seq(3.5, 9.5, 2),
           ymin = 0, ymax = 1, fill = "#b3b3b3", alpha = 0.1)+
  geom_linerange(aes(x = phase,
                     ymin = 0, ymax = Accuracy, group = Model),
                 position = position_dodge(0.6),
                 linetype = "dashed")+
  geom_point(aes(x = phase, y = Accuracy,
                 fill = Model),
             position = position_dodge(0.6),
             size = 3, shape = 21)+
  geom_text(aes(x = phase,
                y = Accuracy+0.05, group = Model,
                label = paste0(round(Accuracy, 2)*100, "%")),
            position = position_dodge(0.6),
            size = 1.5)+
  scale_x_continuous(breaks = 3:10, expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  scale_fill_brewer(name = "", palette = "Set2")+
  xlab("Number of Phases")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        # legend.key.height = unit(1, "mm"),
        # legend.key.width = unit(1, "mm"),
        legend.text = element_text(size = 6, face = "bold.italic"),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, face = "bold.italic"),
        axis.title = element_text(size = 8, face = "bold.italic"))
p_A

# C:
load("Rdata/ensemble_model_9.Rdata")

group_test_data <- str_split(rownames(data_test_meta), "_", simplify = T)[,1]
table(group_test_data)
dataset_names <- unique(group_test_data)

pred_all <- predict(meta_model, newdata = meta_test_matrix)
table(pred_all, labels_test_meta)

acc_datasets <- data.frame(Datasets = dataset_names,
                           Accuracy = NA)
for (i in 1:length(dataset_names)) {
  dataset_ind <- which(group_test_data == dataset_names[i])
  confusion_mat <- as.matrix(table(pred_all[dataset_ind], labels_test_meta[dataset_ind]))
  acc <- round(sum(diag(confusion_mat))/sum(confusion_mat), 4)

  acc_datasets$Accuracy[i] <- acc
}

library(RColorBrewer)

p_C <- ggplot(acc_datasets)+
  geom_linerange(aes(x = Datasets, Accuracy,
                     ymin = 0, ymax = Accuracy), linetype = "dashed")+
  geom_point(aes(Datasets, Accuracy, fill = Datasets),
             size = 6, shape = 21)+
  geom_text(aes(Datasets, Accuracy+0.08,
                label = paste0(round(Accuracy, 4)*100, "%")),
            size = 2)+
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(12))+
  ylim(0, 1)+
  xlab("Names of DataSets")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, face = "bold.italic"),
        axis.title = element_text(size = 8, face = "bold.italic"),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold.italic",
                                  size = 6))

p_C


plot_grid(p_A, NULL, p_C, ncol = 1, labels = c("A", "B", "C"))

