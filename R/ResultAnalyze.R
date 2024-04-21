#' Title Method evaluation: Running time
#'
#' @param func The function that wants to evaluate running time.
#' @param ... The param list of function.
#' @return The running time score of the method.
#' @export
#' @examples
#' data(KICH)
#' data(KEGGgenesetsID)
#' KICHfilter <- NormalizeExp(KICH)
#' RunningTimeScore <- GetRunningTimeScore(UseGSVA, KICHfilter, KEGGgenesetsID)
# GetRunningTimeScore:main function
GetRunningTimeScore <- function(func = NA, ...){
  start_time <- as.numeric(Sys.time())
  result = func(...)
  end_time <- as.numeric(Sys.time())
  return (1.0 / (end_time - start_time))
}

#' Title Method evaluation: Reproducibility score
#'
#' @param exp_feature The gene expression matrix.
#' @param iES_mat The iES_mat obtained by pathway activation evaluation algorithms.
#' @return The reproducibility score.
#' @export
#' @examples
#' data(KICH)
#' data(KEGGgenesetsID)
#' KICHfilter <- NormalizeExp(KICH)
#' pathway_ret <- UseGSVA(KICHfilter, KEGGgenesetsID)
#' ReproducibilityScore <- GetReproducibilityScore(KICHfilter, pathway_ret)
# GetReproducibilityScore: main function
GetReproducibilityScore <- function(exp, iES_mat){
  sample_num <- ncol(exp)
  bp_fun_exp <- function(idx1){
    ret <- c()
    for (idx2 in (idx1+1):sample_num){
      data1 <- exp[,idx1]
      data2 <- exp[,idx2]
      ret <- c(ret, lsa::cosine(data1, data2))
    }
    return (ret)
  }

  bp_fun_path <- function(idx1){
    ret <- c()
    for (idx2 in (idx1+1):sample_num){
      data1 <- iES_mat[,idx1]
      data2 <- iES_mat[,idx2]
      ret <- c(ret, lsa::cosine(data1, data2))
    }
    return (ret)
  }

  results_exp <- BiocParallel::bplapply(1:(sample_num-1), bp_fun_exp)
  results_path <- BiocParallel::bplapply(1:(sample_num-1), bp_fun_path)

  mse_var <- 0
  for (i in 1:(sample_num-1)){
    for (j in 1:(sample_num-i)){
      mse_var <- mse_var + (results_exp[[i]][j] - results_path[[i]][j]) ** 2
    }
  }
  mse_var <- mse_var / ((sample_num*(sample_num-1))/2)
  ret <- 1 / mse_var
  return(ret)
}

#' Title  Method evaluation: ROC AUC
#'
#' @param sample_feature A data frame. The row subscript is the index of the sample. The column subscript is the index of each feature. 
#' @param sample_labels A data frame. The row subscript is the index of the sample. Only one column represents the label.
#' @return ROC AUC score.
#' @export
#' @examples
#' data(KICH)
#' data(KICHcli)
#' data(KEGGgenesetsID)
#' KICHfilter <- NormalizeExp(KICH)
#' pathway_ret <- UseCORGs(KICHfilter, KEGGgenesetsID)
#' ROCAUCscore <- GetROCAUCscore(as.data.frame(t(pathway_ret)), as.data.frame(t(KICHcli[3,])))
# GetROCAUCscore: main function
GetROCAUCscore <- function(sample_feature, sample_labels){
  sample_feature <- as.data.frame(sample_feature)
  sample_labels <- as.data.frame(sample_labels)
  
  cv_num <- 5
  epoch_num <- 50
  colnames(sample_feature) <- paste0("pathway", 1:ncol(sample_feature))
  colnames(sample_labels)[1] <- "labels"
  all_data <- cbind(sample_labels, sample_feature)
  rownames(all_data) <- paste0("sample", 1:nrow(all_data))
  AUC_values <- c()
  for (epoch_cnt in 1:epoch_num){
    p_df <- all_data[all_data$labels == 1, ]
    n_df <- all_data[all_data$labels == 0, ]
    p_folds <- caret::createFolds(1:nrow(p_df), k = cv_num)
    n_folds <- caret::createFolds(1:nrow(n_df), k = cv_num)
    cv_values <- c()
    for (cv_cnt in 1:cv_num) {
      train_p <- p_df[-p_folds[[cv_cnt]], ]
      train_n <- n_df[-n_folds[[cv_cnt]], ]
      test_p <- p_df[p_folds[[cv_cnt]], ]
      test_n <- n_df[n_folds[[cv_cnt]], ]
      train <- rbind(train_p, train_n)
      test <- rbind(test_p, test_n)
      
      model <- randomForest::randomForest(factor(labels) ~ ., data = train, levels = c(unique(sample_labels)))
      pred_probs <- predict(model, newdata = test, type = "prob")
      roc_curve <- pROC::roc((test$labels), pred_probs[, 2])
      cv_values <- c(cv_values, pROC::auc(roc_curve))
    }
    AUC_values <- c(AUC_values, mean(cv_values))
  }
  return (mean(AUC_values))
}