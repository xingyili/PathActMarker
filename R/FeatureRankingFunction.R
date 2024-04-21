


#' Title feature importance ranking obtained by machine learning.
#'
#' @description This function uses machine learning algorithm to get the feature importance ranking based on caret package.
#' @param iES_mat the iES_mat obtained by pathway activation evaluation algorithms
#' @param label normal(0) and tumor (1) in patient samples.
#' @param type Machine learning algorithms, including Random Forest, Linear Regression, Logistic Regression, k-Nearest Neighbors, etc., default random forest.
#' @importFrom  caret trainControl
#' @importFrom caret train
#' @importFrom  caret varImp
#' @return The feature importance ranking
#' @export
#'
#' @examples
#' data(GSVAResults)
#' data(KICHcli)
#' PathwayImportance = FeatureImportance(GSVAResults, label = KICHcli[3,] )
FeatureImportance <-function(iES_mat,label,type = "rf"){
  iES_mat <- t(rbind(iES_mat, label))

  #using the  5-fold cross verification
  control <- caret::trainControl(method="cv", number=5)
  # Training model
  model <- caret::train(as.factor(Label)~., data=iES_mat, method=type, trControl=control)
  # Assess the variable importance
  importance <- caret::varImp(model, scale=FALSE)


  importance = importance$importance
  importance$Variable  <- rownames(importance)
  importance <- importance[order(-importance$Overall),]
  importance_sorted <- as.data.frame(importance[, 1, drop = FALSE])

  return(importance_sorted)

}




#' Title Significance ranking of survival differences between disturbed and normal-like  samples
#'
#'@description This function uses survival analysis to statistically analyze the survival status of disturbed and normal-like  samples in different pathways
#' @param iES_mat the iES_mat obtained by pathway activation evaluation algorithms
#' @param cli the survival information for patient samples.
#' @param samplenum minimum number of normal-like  and perturbed samples
#' @import mclust
#' @import survival
#' @importFrom survminer surv_pvalue
#' @return The survival difference(p-value) between normal-likeand and disturbed samples
#' @export
#'
#' @examples
#' data(GSVAResults)
#' data(KICHcli)
#' PathwayImportance = FeatureImportanceKM(GSVAResults, KICHcli )
FeatureImportanceKM = function(iES_mat, cli, samplenum = 10) {
  cli = as.data.frame(t(cli))

  npaths = nrow(iES_mat)
  path_names = rownames(iES_mat)

  indVec = cli$Label
  inds1 = which(indVec == 0)
  inds2 = which(indVec == 1)
  iES_res = list(normal = iES_mat[, inds1], tumor = iES_mat[, inds2])

  norm_Y = iES_mat[, inds1]
  tumor_Y = iES_mat[, inds2]
  tumor_pat_names = colnames(tumor_Y)
  tumor_com = intersect(tumor_pat_names, row.names(cli))
  tumor_Y = tumor_Y[, which(tumor_pat_names %in% tumor_com)]


  result_df <- data.frame( col_name = character(),  value = numeric(), stringsAsFactors = FALSE )

  for (i in 1:(nrow(tumor_Y))) {
    norm_vec = norm_Y[i,]
    tumor_vec = tumor_Y[i,]
    tmp_m = mclust::Mclust(norm_vec, parameter = TRUE, modelNames = "V")
    id = which.max(tmp_m$parameters$pro)
    tmp_mean = tmp_m$parameters$mean[id]
    tmp_sd = sqrt(tmp_m$parameters$variance$sigmasq[id])

    # UP or DOWN Regulate
    if (tmp_mean < mean(tumor_vec)) {
      thre =  tmp_mean + 2 * tmp_sd
      perturb = names(tumor_vec)[which(tumor_vec >= thre)]
    } else{
      thre =  tmp_mean - 2 * tmp_sd
      perturb = names(tumor_vec)[which(tumor_vec < thre)]
    }
    normlike = setdiff(tumor_com, perturb)
    nlen1 = length(perturb)
    nlen2 = length(normlike)

    if (nlen1 >= samplenum & nlen2 >= samplenum) {
      tmp_cli = rbind(cli[which(row.names(cli) %in% perturb), ],
                      cli[which(row.names(cli) %in% normlike), ])
      tmp_cli$class = c(rep("perturb", nlen1), rep("normlike", nlen2))

      sfit <- survival::survfit(survival::Surv(OS.time, OS) ~ class, data = tmp_cli)
      p_value = survminer::surv_pvalue(sfit, data = tmp_cli)

      row_data <- data.frame(col_name = row.names(tumor_Y)[i], value = p_value$pval.txt)
      result_df <- rbind(result_df, row_data)
    }

  }

  result_df <- result_df[order(result_df$value),]
  row.names(result_df) <- result_df[, 1]
  result_df <- as.data.frame(result_df[, 2, drop = FALSE])

  return(result_df)


}




























