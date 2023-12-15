# 安装和加载必要的包
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

# 载入数据：
load("Rdata/exp_labels.Rdata")
use_python("D:\\Program Files\\Anaconda3", required = T)

# 划分元训练集（70%）和元测试集（30%）
set.seed(123)

index_train_meta <- createDataPartition(y_all, p = 0.7, list = FALSE)
# 确保训练集样本数是5的倍数，并随机打乱：
length(index_train_meta)
index_train_meta <- index_train_meta[sample(1:(length(index_train_meta)-2))]

data_train_meta <- data[index_train_meta, ]
labels_train_meta <- y_all[index_train_meta]
data_test_meta <- data[-index_train_meta, ]
labels_test_meta <- y_all[-index_train_meta]

# 定义基模型的交叉验证函数
cv_folds <- 5
cv <- createFolds(labels_train_meta, k = cv_folds, list = TRUE)

# 初始化基模型的预测结果矩阵
meta_train_mat <- list()
base_model_all <- list()

# 遍历每个基模型
for (i in 1:cv_folds) {
  base_model <- list()

  # 划分基训练集和基验证集
  index_base_train <- unlist(cv[-i])
  index_base_val <- cv[[i]]
  data_base_train <- data_train_meta[index_base_train, ]
  labels_base_train <- labels_train_meta[index_base_train]
  data_base_val <- data_train_meta[index_base_val, ]
  labels_base_val <- labels_train_meta[index_base_val]

  # 初始化基模型的预测结果矩阵
  meta_train_matrix <- matrix(NA, nrow = length(index_base_val), ncol = 0)

  ############ 训练XGBoost模型 ###############
  xgb_model <- xgboost(data = data_base_train,
                       label = as.integer(labels_base_train)-1,
                       nrounds = 500,
                       max_depth = 4,
                       eta = 0.1,
                       gamma = 0,
                       colsample_bytree = 0.8,
                       min_child_weight = 1,
                       subsample = 0.7,
                       objective = "multi:softprob",
                       num_class = 9)
  base_model[[1]] <- xgb_model
  xgb_predictions <- as.numeric(matrix(predict(xgb_model, data_base_val),
                                       ncol = 9, byrow = TRUE) %>% k_argmax())
  # 看看准确率情况：
  confusion_mat <- table(xgb_predictions, labels_base_val)
  message("XGBoost: Validation ", i,
          "\nAccuracy: ",
          round(sum(diag(confusion_mat))/sum(confusion_mat), 3))

  # 将验证集的结果存入矩阵：
  meta_train_matrix <- cbind(meta_train_matrix, as.vector(xgb_predictions))

  ############ 训练SVM模型 ############
  svm_model <- svm(labels_base_train ~ ., data = data_base_train)
  base_model[[2]] <- svm_model
  svm_predictions <- predict(svm_model, newdata = data_base_val)
  # 看看准确率情况：
  confusion_mat <- table(svm_predictions, labels_base_val)
  message("SVM: Validation ", i,
          "\nAccuracy: ",
          round(sum(diag(confusion_mat))/sum(confusion_mat), 3))
  # 将验证集的结果存入矩阵：
  meta_train_matrix <- cbind(meta_train_matrix, as.vector(svm_predictions))


  ############ 训练随机森林模型 ############
  rf_model <- randomForest(labels_base_train ~ ., data = data_base_train, ntree = 100)
  base_model[[3]] <- rf_model
  rf_predictions <- predict(rf_model, newdata = data_base_val, type = "response")
  # 看看准确率情况：
  confusion_mat <- table(rf_predictions, labels_base_val)
  message("RF: Validation ", i,
          "\nAccuracy: ",
          round(sum(diag(confusion_mat))/sum(confusion_mat), 3))
  # 将验证集的结果存入矩阵：
  meta_train_matrix <- cbind(meta_train_matrix, as.vector(rf_predictions))


  ############ 将验证集的结果存入矩阵 ############
  # 训练神经网络模型
  labels_base_train <- to_categorical(labels_base_train, 9)

  # 创建一个Sequential模型
  nn_model <- keras_model_sequential()

  # 添加神经网络层
  # 假设你有三个基因的表达量数据，每个基因作为一个特征输入
  nn_model %>%
    layer_dense(units = 64, activation = "sigmoid", input_shape = c(152)) %>%
    layer_dense(units = 32, activation = "sigmoid") %>%
    layer_dense(units = 16, activation = "sigmoid") %>%
    layer_dense(units = 9, activation = "softmax")  # 输出层，使用 softmax 以获得分类概率

  # 编译模型
  nn_model %>% compile(
    loss = "categorical_crossentropy",
    optimizer = optimizer_adam(learning_rate = 0.001),
    metrics = c("accuracy")
  )

  # 训练模型
  nn_model %>% fit(
    x = data_base_train,
    y = labels_base_train,
    epochs = 50,
    batch_size = 32,
    validation_split = 0
  )

  save_model_tf(nn_model, paste0("Rdata/nn_models/model_", i))
  # base_model[[4]] <- nn_model
  nn_predictions <- as.numeric(predict(nn_model, data_base_val) %>% k_argmax())
  # 看看准确率情况：
  confusion_mat <- table(nn_predictions, labels_base_val)
  message("Neural Network: Validation ", i,
          "\nAccuracy: ",
          round(sum(diag(confusion_mat))/sum(confusion_mat), 3))

  # 将验证集的结果存入矩阵：
  meta_train_matrix <- cbind(meta_train_matrix, as.vector(nn_predictions))

  meta_train_mat[[i]] <- meta_train_matrix
  base_model_all[[i]] <- base_model
}

# 合并5次交叉验证的结果作为元模型训练集：
meta_train_matrix <- Reduce(rbind, meta_train_mat)
labels_train_meta2 <- labels_train_meta[unlist(cv)]
colnames(meta_train_matrix) <- c("xgboost", "svm", "rf", "nn")

############### 定义元模型（这里使用随机森林）##########
meta_model <- randomForest(labels_train_meta2 ~ ., data = meta_train_matrix)

# 元模型在训练集上的表现：
stack_predictions <- predict(meta_model, newdata = meta_train_matrix)
confusionMatrix(stack_predictions, labels_train_meta2)

# 元模型在测试集上的预测
# 通过get_meta_test_mat函数构建测试集：
source("get_meta_test_mat.R")
meta_test_matrix <- get_meta_test_mat(base_model_all, data_test_meta)

final_predictions <- predict(meta_model,
                             newdata = meta_test_matrix)

# 评估最终模型的性能
confusionMatrix(final_predictions, labels_test_meta)
# 集成后模型准确率相对最高的神经网络提高了3%，达到了81%左右。

save(base_model_all, meta_model, data_test_meta,
     file = "../细胞阻滞探究/Rdata/ensemble_model.Rdata")
