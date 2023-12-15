get_meta_test_mat <- function(base_model_all, data_test_meta,
                              nn_model_dir = "Rdata/nn_models",
                              label_test_meta = NULL,
                              phase_num = 9){
  ########### 平均化5个模型的测试结果作为测试集 ############
  pred_test_all <- list()
  for (i in 1:5) {
    xgb_pred_test <- as.numeric(matrix(predict(base_model_all[[i]][[1]], data_test_meta),
                                       ncol = phase_num, byrow = TRUE) %>% k_argmax())
    message("xgb ", i, " done!")
    svm_pred_test <- predict(base_model_all[[i]][[2]], data_test_meta)
    message("svm ", i, " done!")
    rf_pred_test <- predict(base_model_all[[i]][[3]], data_test_meta)
    message("rf ", i, " done!")
    nn_model <- load_model_tf(paste0(nn_model_dir, "/model_", i))
    nn_pred_test <- as.numeric(predict(nn_model, data_test_meta) %>% k_argmax())
    message("nn ", i, " done!")

    if (!is.null(label_test_meta)) {
      ###### 只在测试集有标签的时候运行 ------------
      # XGB准确率：
      confusion_mat <- table(xgb_pred_test, labels_test_meta)
      message("XGBoost: Test ", i,
              "\nAccuracy: ",
              round(sum(diag(confusion_mat))/sum(confusion_mat), 3))
      # SVM准确率：
      confusion_mat <- table(svm_pred_test, labels_test_meta)
      message("\nSVM: Test ", i,
              "\nAccuracy: ",
              round(sum(diag(confusion_mat))/sum(confusion_mat), 3))
      # RF准确率：
      confusion_mat <- table(rf_pred_test, labels_test_meta)
      message("\nRF: Test ", i,
              "\nAccuracy: ",
              round(sum(diag(confusion_mat))/sum(confusion_mat), 3))
      # NN准确率：
      confusion_mat <- table(nn_pred_test, labels_test_meta)
      message("\nNeural Network: Test ", i,
              "\nAccuracy: ",
              round(sum(diag(confusion_mat))/sum(confusion_mat), 3))
    }

    pred_test <- cbind(xgb_pred_test, as.vector(svm_pred_test),
                       as.vector(rf_pred_test), nn_pred_test)
    pred_test_all[[i]] <- pred_test
  }

  # 获取向量中出现最多的元素：
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  meta_test_matrix <- matrix(NA, nrow = nrow(data_test_meta), ncol = 0)

  for (i in 1:4) {
    pred_model5 <- Reduce(cbind, lapply(pred_test_all, function(x) x[,i]))
    pred_mean <- apply(pred_model5, 1, getmode)
    meta_test_matrix <- cbind(meta_test_matrix, pred_mean)
  }

  colnames(meta_test_matrix) <- c("xgboost", "svm", "rf", "nn")

  return(meta_test_matrix)
}
