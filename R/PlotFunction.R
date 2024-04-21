#' Title Plot Radar
#'
#' @param performanceScores a data frame, rows are all methods, cols are all scores.
#' @return The radar plot.
#' @export
#'
#' @examples
#' data(KICH)
#' data(KEGGgenesetsID)
#' data(KICHcli)
#' KICHfilter <- NormalizeExp(KICH)
#' GSVARunningTimeScore <- GetRunningTimeScore(UseGSVA, KICHfilter, KEGGgenesetsID)
#' GSVApathway_ret <- UseGSVA(KICHfilter, KEGGgenesetsID)
#' GSVAReproducibilityScore <- GetReproducibilityScore(KICHfilter, GSVApathway_ret)
#' GSVAROCAUCscore <- GetROCAUCscore(as.data.frame(t(GSVApathway_ret)), as.data.frame(t(KICHcli[3,])))
#' CORGsRunningTimeScore <- GetRunningTimeScore(UseCORGs, KICHfilter, KEGGgenesetsID)
#' CORGspathway_ret <- UseCORGs(KICHfilter, KEGGgenesetsID)
#' CORGsReproducibilityScore <- GetReproducibilityScore(KICHfilter, CORGspathway_ret)
#' CORGsROCAUCscore <- GetROCAUCscore(as.data.frame(t(CORGspathway_ret)), as.data.frame(t(KICHcli[3,])))
#' ssGSEARunningTimeScore <- GetRunningTimeScore(UsessGSEA, KICHfilter, KEGGgenesetsID)
#' ssGSEApathway_ret <- UsessGSEA(KICHfilter, KEGGgenesetsID)
#' ssGSEAReproducibilityScore <- GetReproducibilityScore(KICHfilter, ssGSEApathway_ret)
#' ssGSEAROCAUCscore <- GetROCAUCscore(as.data.frame(t(ssGSEApathway_ret)), as.data.frame(t(KICHcli[3,])))
#' PlotRadar(data.frame(t(data.frame(row.names = c("RunningTimeScore", "ReproducibilityScore", "ROCAUCscore"), 
#'                                  GSVA = c(GSVARunningTimeScore, GSVAReproducibilityScore, GSVAROCAUCscore),
#'                                  CORGs = c(CORGsRunningTimeScore, CORGsReproducibilityScore, CORGsROCAUCscore),
#'                                  ssGSEA = c(ssGSEARunningTimeScore, ssGSEAReproducibilityScore, ssGSEAROCAUCscore)))))
PlotRadar = function(performanceScores){
  colors <- tail(RColorBrewer::brewer.pal(nrow(performanceScores) + 2, "Greys"), nrow(performanceScores))
  # colors <- grDevices::colorRampPalette(c("darkblue", "blue"))(nrow(performanceScores))
  colors_border <- colors
  colors_in <- scales::alpha(colors,0.3)
  
  fmsb::radarchart(performanceScores  , axistype=0 , maxmin=F,
              #polygon
              pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
              #the grid
              cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
              #labels
              vlcex=0.8
  )
  legend(x=0.7, y=1, legend = rownames(performanceScores), 
         bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
}


#' Title Plot Heatmap
#'
#' @param iES_mat the iES_mat obtained by pathway activation evaluation algorithms.
#' @param label normal(0) and tumor (1) in patient samples.
#' @param condition boolean value indicating whether samples need to be grouped and colored according to different conditions.
#' @param title a string representing the title of the heat map.
#' @import pheatmap
#' @return p
#' @export
#'
#' @examples
#' data(GSVAResults)
#' data(KICHcli)
#' PlotHeatmap(GSVAResults, label = KICHcli[3,])
PlotHeatmap = function(iES_mat, label, condition = TRUE ,title = NA){

  myColor = colorRampPalette(c("Darkblue", "white","red"))(100)
  iES_matorder = iES_mat[, order(as.numeric(label))]

  if (condition == TRUE){
    normalvalue <- sum(label == 0)
    diseasevalue <- sum(label == 1)
    annotation_col = data.frame(condition = c(rep("Normal Samples",normalvalue),rep("Diease Samples",diseasevalue)))  #增加Time，CellType分组信息
    rownames(annotation_col) =  colnames(iES_matorder)

    p =  pheatmap::pheatmap(iES_matorder, annotation_col = annotation_col,show_colnames = FALSE,cluster_rows = F,cluster_cols = F,
                  color=myColor, main = title, angle_col = 45, treeheight_col = 0,
                  border_color = NA)

  }else{
    p = pheatmap::pheatmap(iES_matorder,show_colnames = FALSE,cluster_rows = F,cluster_cols = F,
                 color=myColor, main = title, angle_col = 45, treeheight_col = 0,
                 border_color = NA)
  }
  return(p)
}




#' Title  Plot Waterfall
#'
#' @param iES_mat the iES_mat obtained by pathway activation evaluation algorithms
#' @param path_name  the pathway name
#' @param label normal(0) and tumor (1) in patient samples.
#' @param title  boolean true or false for including the title (path_name) in the Waterfall plot.
#'
#' @import ggplot2
#' @import ggpubr
#' @return p
#' @export
#'
#' @examples
#' data(GSVAResults)
#' data(KICHcli)
#' PlotWaterfall(GSVAResults, path_name = "HALLMARK_PI3K_AKT_MTOR_SIGNALING" , label = KICHcli[3,] )
PlotWaterfall = function(iES_mat, path_name,label, title = TRUE){
  group_colors = c(tumor = "Brown", normal = "#56B4E9")
  sort1 = which(label==0)
  sort2 = which(label==1)
  iES_res = list(normal = iES_mat[, sort1], tumor = iES_mat[, sort2])

  mat_sort = which(rownames(iES_mat) == path_name)
  normal = iES_res[[1]][mat_sort,]; tumor = iES_res[[2]][mat_sort, ]

  n_gap = round((length(tumor) + length(normal)) * 0.01)
  tmp_iES_mat = data.frame(value = c(tumor[order(tumor)], rep(0, n_gap), normal[order(normal)]),
                           type = c(rep("tumor", length(tumor)), rep("tumor", n_gap), rep("normal", length(normal))),
                           fill = c(rep("fill", length(tumor)), rep(NA, n_gap), rep("fill", length(normal))))
  nrow_iES_mat = nrow(tmp_iES_mat)
  if (title ==TRUE){
    p = ggplot(tmp_iES_mat, aes(x = seq_len(nrow_iES_mat), y = value)) +
      geom_area(aes(fill = type)) +
      theme(legend.position="top", legend.direction="horizontal", panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      labs(x="Samples", y="Enrichment Score") +
      scale_fill_manual(values=group_colors) +
      ggtitle(path_name)
  }else{
    p = ggplot(tmp_iES_mat, aes(x = seq_len(nrow_iES_mat), y = value)) +
      geom_area(aes(fill = type)) +
      theme(legend.position="top", legend.direction="horizontal", panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      labs(x="Samples", y="Enrichment Score") +
      scale_fill_manual(values=group_colors)
  }
  return(p)
}




#' Title Plot densityfall
#'
#' @param iES_mat the iES_mat obtained by pathway activation evaluation algorithms
#' @param path_name the pathway name
#' @param label normal(0) and tumor (1) in patient samples.
#' @param title boolean true or false for including the title (path_name) in the densityfall plot.
#' @import ggpubr
#' @import mclust
#' @importFrom survminer ggsurvplot
#' @return p
#' @export
#'
#' @examples
#' data(GSVAResults)
#' data(KICHcli)
#' PlotDensityfall(GSVAResults, path_name = "HALLMARK_PI3K_AKT_MTOR_SIGNALING" , label = KICHcli[3,] )
PlotDensityfall = function(iES_mat, path_name,label, title = TRUE){
  sort1 = which(label==0)
  sort2 = which(label==1)
  iES_res = list(normal = iES_mat[, sort1], tumor = iES_mat[, sort2])

  mat_sort = which(rownames(iES_mat) == path_name)
  normal = iES_res[[1]][mat_sort,]; tumor = iES_res[[2]][mat_sort, ]
  tmp_m = Mclust(normal, parameter = TRUE, modelNames = "V")
  id = which.max(tmp_m$parameters$pro)
  tmp_mean = tmp_m$parameters$mean[id]
  tmp_iES_mat = data.frame(value = c(normal, tumor),
                           type = c(rep("normal", length(normal)), rep("tumor", length(tumor))))
  if (title ==TRUE){
    p = ggdensity(tmp_iES_mat, x = "value",rug = TRUE,
                  color = "type", fill = "type",
                  palette = c("#56B4E9", "Brown"),
                  main = path_name, legend.title = "") +
      geom_vline(xintercept= c(tmp_mean, mean(tumor)), linetype="dashed", color = c("#56B4E9", "Brown"))
  }else{
    p = ggdensity(tmp_iES_mat, x = "value", rug = TRUE,
                  color = "type", fill = "type",
                  palette = c("#56B4E9", "Brown"),
                  legend.title = "")+
      geom_vline(xintercept= c(tmp_mean, mean(tumor)), linetype="dashed", color = c("#56B4E9", "Brown"))
  }
  return(p)

}



#' Title Plot Kaplan-Meier(KM) survival curve
#'
#' @param iES_mat the the iES_mat obtained by pathway activation evaluation algorithms.
#' @param cli the survival information for patient samples.
#' @param path_name the pathway name.
#' @param sample_num the minimum number of patients for survival analysis.
#' @param title boolean true or false for including the title (path_name) in the KM survival plot.
#' @importFrom survminer ggsurvplot
#' @import survival
#' @import mclust
#' @import ggpubr
#' @return p
#' @export
#'
#' @examples
#' data(GSVAResults)
#' data(KICHcli)
#' PlotiESSurv(GSVAResults, cli = KICHcli, path_name = "HALLMARK_PI3K_AKT_MTOR_SIGNALING")
PlotiESSurv = function(iES_mat,cli, path_name, sample_num = 20, title = TRUE) {
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

  norm_vec = norm_Y[path_name,]
  tumor_vec = tumor_Y[path_name,]
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

  if (nlen1 >= sample_num & nlen2 >= sample_num) {
    tmp_cli = rbind(cli[which(row.names(cli) %in% perturb), ],
                    cli[which(row.names(cli) %in% normlike), ])
    tmp_cli$class = c(rep("perturb", nlen1), rep("normlike", nlen2))

    tmp_cox = survival::coxph(survival::Surv(OS.time, OS) ~ class,
                              data = tmp_cli,
                              method = "breslow")
    nperturb = length(perturb)
    sfit <-
      survival::survfit(survival::Surv(OS.time, OS) ~ class,
                        data = tmp_cli,
                        conf.type = "log-log")

    if (title == TRUE) {
      p = survminer::ggsurvplot(
        sfit,
        conf.int = TRUE,
        pval = TRUE,
        risk.table = FALSE,
        legend.labs = c("perturbed", "normal-like"),
        legend.title = "",
        palette = c("darkred", "darksalmon"),
        linetype = "strata",
        title = path_name,
        data = tmp_cli
      )

    } else{
      p = survminer::ggsurvplot(
        sfit,
        conf.int = TRUE,
        pval = TRUE,
        risk.table = FALSE,
        legend.labs = c("perturbed", "normal-like"),
        legend.title = "",
        palette = c("darkred", "darksalmon"),
        linetype = "strata",
        data = tmp_cli
      )
    }
  } else{
    p = NULL
  }
  return(p)
}






#' Title  Obtain adjusted rand index (ARI) score
#'
#' @param iES_mat the iES_mat obtained by pathway activation evaluation algorithms.
#' @param label normal(0) and tumor (1) in patient samples.
#' @param samplenum the number of tumor and normal samples
#' @param n_iter number of repeated experiments
#' @import limma
#' @import cluster
#' @import mclust
#' @importFrom dendextend cutree
#' @return ARI_vec
#' @export
#'
#' @examples
#' data(GSVAResults)
#' data(KICHcli)
#' GSVAARIscore = GetARIscore(GSVAResults,label = KICHcli[3,] )
GetARIscore = function(iES_mat, label, samplenum = 20,n_iter = 50) {


  ARI_vec <- numeric(n_iter)
  for (i in 1:n_iter) {

    normal_sample <- sample(which(label[1,] == 0), samplenum, replace = F)
    tumor_sample <- sample(which(label[1,] == 1), samplenum, replace = F)
    selected_sample <- c(normal_sample, tumor_sample)
    iES_mat_selected <- as.matrix(iES_mat[,selected_sample])


    sample_types <- factor(c(rep("normal", samplenum), rep("tumor", samplenum)))
    design <- model.matrix(~ sample_types)
    fit <- limma::lmFit(iES_mat_selected, design)
    fit <- limma::eBayes(fit)
    results <- limma::topTable(fit, coef = 2, n = Inf, sort.by = "p", adjust.method = "BH")

    # select the top 10 gene sets according to the adjusted p values and perform the hierarchical clustering.
    top_gene_sets <- head(order(fit$p.value[,2]), 10)
    top_gene_set_es <- iES_mat_selected[top_gene_sets, ]
    hc <- hclust(dist(t(top_gene_set_es)), method = "ward.D2")
    dend <- as.dendrogram(hc)

    # bipartition the hierarchical tree into two classes and compare the clustering results with sample labels using ARI.
    clusters <- dendextend::cutree(dend, k = 2)
    clusters<- ifelse(clusters == 1, 0, ifelse(clusters == 2, 1, clusters))
    ARI <- mclust::adjustedRandIndex(as.numeric(label[1, selected_sample]), as.numeric(clusters))
    ARI_vec[i] <- ARI

  }

  # return(mean(ARI_vec))
  return (ARI_vec)

}




#' Title Plot ggviolin
#'
#' @param ARIscores a two-column matrix, the first column(methods) and the second column(ARI value)
#' @param compare Boolean value indicating whether statistical tests are performed for comparison between groups.
#' @import  ggpubr
#' @return p
#' @export
#'
#' @examples
#' data(ARIscores)
#' plotggviolin(ARIscores)
plotggviolin = function(ARIscores, compare = F) {
  x = colnames(ARIscores)[1]
  y = colnames(ARIscores)[2]
  colnames(ARIscores) = c("Methods", "ARIs")


  group = as.character(unique(factor(ARIscores$Methods)))
  ARIscores$Methods = factor(ARIscores$Methods, levels = group)
  comp = combn(group, 2)
  my_comparisons = list()
  for (i in 1:ncol(comp)) {

    my_comparisons[[i]] <- comp[, i]
    if (i>=4){
      break
    }
  }

  if (compare == TRUE) {
    p =  ggpubr::ggviolin( ARIscores, x = "Methods", y = "ARIs", fill = "Methods",
                           xlab = F, ylab = y,  legend.title = x, add = "boxplot",  add.params = list(fill = "white"),
                           scale = "width",  width = 0.6  ) + ggpubr::stat_compare_means(comparisons = my_comparisons) + coord_cartesian(ylim = c(0, 1.5)) +
      scale_y_continuous(breaks = c(0, 0.5, 1))+  theme(legend.key.size = unit(0.8, "lines"))

  } else{
    p =  ggpubr::ggviolin( ARIscores, x = "Methods", y = "ARIs", fill = "Methods",
                           xlab = F, ylab = y,  legend.title = x, add = "boxplot",  add.params = list(fill = "white"),
                           scale = "width",  width = 0.6  )
  }

  return(p)

}