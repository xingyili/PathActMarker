

#' Title Filter out genes with constant expression values
#' @description This is an R function that filters out genes with constant expression values from an expression matrix.
#' Genes with standard deviation less than a threshold value of 1e-10 are considered to have constant expression values
#' and are removed from the expression matrix.
#' @param expr Gene expression matrix
#' @keywords internal
#' @importFrom matrixStats rowSds
#' @return expr
#'
.filterFeatures <- function(expr) {
  expr <- as.matrix(expr)
  sdGenes <- matrixStats::rowSds(expr)

  sdGenes[sdGenes < 1e-10] <- 0
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(
      sum(sdGenes == 0 | is.na(sdGenes)),
      " genes with constant expression values throuhgout the samples."
    )

    expr <- expr[sdGenes > 0 & !is.na(sdGenes),]

  }

  if (nrow(expr) < 2)
    stop("Less than two genes in the input assay object\n")

  if (is.null(rownames(expr)))
    stop("The input assay object doesn't have rownames\n")

  return(expr)
}




#' Title Expression matrix standardization
#'@description This function is used to standardize the gene expression matrix, which mainly consists of three steps:
#'Step 1: Filter out the same gene Step 2: Filter out genes with constant expression values Step 3: Normalize the matrix
#' @param exp  Gene expression matrix
#' @param method Z-score, quantile, TPM, TMM, min-max
#' @import edgeR
#' @import limma
#' @return  The normalized matrix
#' @export
#'
#' @examples
#'data(KICH)
#'KICHfilter = NormalizeExp(KICH)
#'
NormalizeExp <- function(exp, method = "none") {
  #step1:Filter out the same gene
  exp <- as.data.frame(exp)
  exp <- exp[!duplicated(exp$GeneId),]

  #step2:Filter out genes with constant expression values
  rownames(exp) <- exp[, 1]
  exp <- exp[, -1]
  expfilter = .filterFeatures(exp)

  #step3:normalization method
  if (method == "TPM") {
    # TPM
    rowsums <- apply(expfilter, 1, sum)
    expnormal <- t(t(expfilter) / rowsums * 1e6)
  }  
  else if (method == "TMM") {
    # TMM
    dge <- edgeR::DGEList(expfilter)
    keep <- rowSums(edgeR::cpm(dge) > 1) >= 2
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dgelist_norm <- edgeR::calcNormFactors(dge, method = 'TMM')
    expnormal = dgelist_norm$counts
    # quantile
  } 
  else if (method == "quantile") {
    expnormal <-limma::normalizeBetweenArrays(expfilter, method = "quantile")
  }
  else if (method == "min-max") {
    # min-max
    expnormal <- expfilter
    for (i in 1:ncol(expfilter)) {
      expnormal[,i] <- (expfilter[,i] - min(expfilter[,i])) / (max(expfilter[,i]) - min(expfilter[,i]))
    }
  }
  else if (method == "z-score"){
    # Z-score
    expnormal <- scale(expfilter)
  }
  else{
    expnormal <- expfilter
  }
  return(expnormal)
}
