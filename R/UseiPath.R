# Pathway activation evaluation algorithm：IPATH




#' remove genes with 0 sd
#'
#' @description This function uses remove non-informative genes.
#' @param exp is the gene expression matrix.
#' @importFrom stats sd
#' @keywords internal
#' @return a processed matrix
remove_data = function(exp){
  rem_ids = which(apply(exp, 1, sd) == 0)
  if (length(rem_ids) == 0){
    return(exp)
  }else{
    return(exp[-rem_ids, ])
  }
}


#' set up for the parallel computing for biocParallel.
#'
#' @description This function sets up the environment for parallel computing.
#' @param nprocess number of processors
#' @param BPPARAM bpparameter from bpparam
#' @keywords internal
#' @return BAPPARAM settings
setUp_BPPARAM = function (nprocess = 0, BPPARAM = NULL)
{
  if (is.null(BPPARAM)) {
    if (nprocess != 0) {
      if (.Platform$OS.type == "windows") {
        result <- SnowParam(workers = nprocess)
      }
      else {
        result <- MulticoreParam(workers = nprocess)
      }
    }
    else {
      result <- BiocParallel::bpparam()
    }
    return(result)
  }
  else {
    return(BPPARAM)
  }
}


#' iES calculation Function
#'
#' This function calculates the iES matrix which is the core of iPath.
#' @importFrom matrixStats rowSds
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit
#' @import BiocParallel
#' @param exp  the gene expression matrix.
#' @param pathwaylist  the pathway genesets.
#' @param BPPARAM parameters from the BiocParallel.
#' @param nPro number of processors (default = 0).
#' @return a matrix with rows corresponding to the pathways and columns corresponding to the patients.
#' @keywords iES statistics calculatioin.
#' @export
#'
#' @examples
#' data(KICH)
#' data(KEGGgenesetsID)
#' KICHfilter <- NormalizeExp(KICH)
#' iPathResults = UseiPath(KICHfilter, KEGGgenesetsID)
UseiPath = function(exp, pathwaylist,  BPPARAM = NULL, nPro = 0){
  exp = as.matrix(exp)
  npats = ncol(exp)
  GSDB_paths = pathwaylist
  GSDB_paths_names = names(pathwaylist)
  ngsets = length(GSDB_paths);
  message("start normalization ...")
  exp = remove_data(exp)
  ngenes = nrow(exp)
  gnames = rownames(exp)
  row_mean = rowMeans(exp)
  row_sd = matrixStats::rowSds(exp)
  tmp_norm = abs((exp-row_mean) / row_sd)

  order_array = apply(tmp_norm, 2, function(x) rev(order(x)))
  order_name = apply(order_array, 2, function(i) gnames[i])
  order_stats = vapply(seq_len(npats), function(i) tmp_norm[, i][order_array[, i]],
                       FUN.VALUE = numeric(ngenes))

  bp_fun = function(i) {
    one_stats = order_stats[, i]
    one_names = order_name[, i]
    names(one_stats) = one_names
    one_pat_vec = vapply(seq_len(ngsets), function(j){
      one_match_pos = na.omit(match(GSDB_paths[[j]], one_names))
      return(caliES2(one_stats, one_match_pos))
    }, FUN.VALUE = numeric(1))
  }

  tmpParam = setUp_BPPARAM(nPro, BPPARAM = BPPARAM)
  tmpParam$progressbar = TRUE
  message("start calculating iES matrix ...")

  pats_iESs = BiocParallel::bplapply(seq_len(npats), FUN = bp_fun, BPPARAM = tmpParam)
  iES_mat = do.call(cbind, pats_iESs)
  dimnames(iES_mat) = list(GSDB_paths_names, colnames(exp))
  return(iES_mat)
}
