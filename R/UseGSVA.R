# Pathway activation evaluation algorithm：GSVA



#' Title  Pathway activation evaluation algorithm：GSVA
#'
#' @param exp the gene expression matrix.
#' @param pathwaylist the pathway genesets.
#' @param kernelMethods character string denoting the kernel to use during the non-parametric estimation of the cumulative distribution function of expression levels across samples
#' when method="gsva". By default, kcdf="Gaussian" which is suitable when input expression values are continuous, such as microarray.
#' When input expression values are integer counts, such as RNA-seq experiments should be set to kcdf="Poisson"
#' @param min.sz minimum size of the resulting gene sets.
#' @param max.sz maximum size of the resulting gene sets.
#' @param parallel.sz number of threads of execution to use when doing the calculations in parallel.
#' @param mx.diff  two approaches to calculate the enrichment statistic (ES) from the KS random walk statistic. mx.diff=FALSE: ES is calculated as the maximum
#'  distance of the random walk from 0. mx.diff=TRUE (default): ES is calculated as the magnitude difference between the largest positive and negative random walk deviations.
#' @param tau exponent defining the weight of the tail in the random walk
#' @param ssgsea.norm Logical, set to TRUE (default) . When ssgsea.norm=FALSE this last normalization     step is skipped.
#' @param verbose Gives information about each calculation step. Default: FALSE.
#' @param BPPARAM An object of class BiocParallelParam specifiying parameters.
#' @param abs.ranking When abs.ranking=FALSE (default) a modified Kuiper statistic is used to calculate enrichment scores.
#' @importFrom BiocParallel SerialParam
#' @return IES
#' @export
#' @examples
#' data(KICH)
#' data(KEGGgenesetsID)
#' KICHfilter <- NormalizeExp(KICH)
#' GSVAResults = UseGSVA(KICHfilter, KEGGgenesetsID)
UseGSVA <-function(exp, pathwaylist,
                   kernelMethods=c("Gaussian", "Poisson", "none"),
                   abs.ranking=FALSE,
                   min.sz=1,
                   max.sz=Inf,
                   parallel.sz=1L,
                   mx.diff=TRUE,
                   tau=1,
                   ssgsea.norm=TRUE,
                   verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)){

  kernelMethods <- match.arg(kernelMethods)

  exp <- .filterFeatures(exp)


  mapped.pathwaylist <- .mapGeneSetsToFeatures(pathwaylist, rownames(exp))

  mapped.pathwaylist <- filterGeneSets(mapped.pathwaylist,
                                       min.sz=max(1, min.sz),
                                       max.sz=max.sz)

  if (!missing(kernelMethods)) {
    if (kernelMethods == "Gaussian") {
      rnaseq <- FALSE
      kernel <- TRUE
    } else if (kernelMethods == "Poisson") {
      rnaseq <- TRUE
      kernel <- TRUE
    } else
      kernel <- FALSE
  }

  IES <- .gsva(exp, mapped.pathwaylist, kernelMethods, rnaseq, abs.ranking,
               parallel.sz, mx.diff, tau, kernel, ssgsea.norm, verbose, BPPARAM)

  return(IES)

}



#' @keywords internal
#' @importFrom BiocParallel SerialParam
#' @importMethodsFrom BiocParallel bpworkers
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel bpworkers
.gsva <- function(exp, pathwaylist,
                  kernelMethods=c("Gaussian", "Poisson", "none"),
                  rnaseq=FALSE,
                  abs.ranking=FALSE,
                  parallel.sz=1L,
                  mx.diff=TRUE,
                  tau=1,
                  kernel=TRUE,
                  ssgsea.norm=TRUE,
                  verbose=TRUE,
                  BPPARAM=SerialParam(progressbar=verbose)) {


  if (length(pathwaylist) == 0)
    stop("The gene set list is empty! Filter may be too stringent.")

  if (any(lengths(pathwaylist) == 1))
    warning("Some gene sets have size one. Consider setting 'min.sz > 1'.")

  parallel.sz <- as.integer(parallel.sz)
  if (parallel.sz < 1L)
    parallel.sz <- 1L

  if (parallel.sz > 1L && class(BPPARAM) == "SerialParam") {
    BPPARAM=MulticoreParam(progressbar=verbose, workers=parallel.sz, tasks=100)
  } else if (parallel.sz == 1L && class(BPPARAM) != "SerialParam") {
    parallel.sz <- bpnworkers(BPPARAM)
  } else if (parallel.sz > 1L && class(BPPARAM) != "SerialParam") {
    bpworkers(BPPARAM) <- parallel.sz
  }

  if (class(BPPARAM) != "SerialParam" && verbose)
    cat(sprintf("Setting parallel calculations through a %s back-end\nwith workers=%d and tasks=100.\n",
                class(BPPARAM), parallel.sz))

  if(verbose)
    cat("Estimating gsva scores for", length(pathwaylist),"gene sets.\n")

  n.samples <- ncol(exp)
  n.genes <- nrow(exp)
  n.gset <- length(pathwaylist)

  es.obs <- matrix(NaN, n.gset, n.samples, dimnames=list(names(pathwaylist),colnames(exp)))
  colnames(es.obs) <- colnames(exp)
  rownames(es.obs) <- names(pathwaylist)

  es.obs <- compute.geneset.es(exp, pathwaylist, 1:n.samples,
                               rnaseq=rnaseq, abs.ranking=abs.ranking,
                               parallel.sz=parallel.sz,
                               mx.diff=mx.diff, tau=tau, kernel=kernel,
                               verbose=verbose, BPPARAM=BPPARAM)

  colnames(es.obs) <- colnames(exp)
  rownames(es.obs) <- names(pathwaylist)

  return(es.obs)
}





















































