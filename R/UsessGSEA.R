# Pathway activation evaluation algorithm：ssGSEA


#' Title  filter out genes with constant expression values
#' @param expr the gene exp
#' @importFrom matrixStats rowSds
.filterFeaturesgsea <- function(expr) {
  expr <- as.matrix(expr)

  sdGenes <- matrixStats::rowSds(expr)
  sdGenes[sdGenes < 1e-10] <- 0
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            " genes with constant expression values throuhgout the samples.")
  }

  if (nrow(expr) < 2)
    stop("Less than two genes in the input assay object\n")

  if(is.null(rownames(expr)))
    stop("The input assay object doesn't have rownames\n")

  expr
}


#' @keywords internal
#' @importMethodsFrom IRanges match
#' @importFrom IRanges CharacterList
.mapGeneSetsToFeatures <- function(gsets, features) {

  gsets2 <- IRanges::CharacterList(gsets)
  mt <- match(gsets2, features)
  mapdgenesets <- as.list(mt[!is.na(mt)])

  if (length(unlist(mapdgenesets, use.names=FALSE)) == 0)
    stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")

  return(mapdgenesets)
}





#' @keywords internal
.sparseToList <-function(dgCMat, MARGIN){
  MARGIN <- as.integer(MARGIN)
  J <- rep(1:ncol(dgCMat), diff(dgCMat@p))
  I <- dgCMat@i + 1
  x <- dgCMat@x
  if (MARGIN == 1L) {
    result <- split(x, I)
    names(result) <- rownames(dgCMat)[as.numeric(names(result))]
  } else if (MARGIN == 2L) {
    result <- split(x, J)
    names(result) <- colnames(dgCMat)[as.numeric(names(result))]
  }
  else {
    warning("invalid MARGIN; return NULL")
    result <- NULL
  }
  return(result)
}



#' @keywords internal
.dgCapply<-function(m,f, MARGIN){
  x <- lapply(.sparseToList(m, MARGIN), f)
  m@x <- unlist(x, use.names=FALSE)
  m
}



#' Title  Pathway activation evaluation algorithm：ssGSEA
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
#'
#' @return IES
#' @export
#' @examples
#' data(KICH)
#' data(KEGGgenesetsID)
#' KICHfilter <- NormalizeExp(KICH)
#' ssGSEAResults = UsessGSEA(KICHfilter, KEGGgenesetsID)
UsessGSEA <-function(exp, pathwaylist,
                     kernelMethods=c("Gaussian", "Poisson", "none"),
                     abs.ranking=FALSE,
                     min.sz=1,
                     max.sz=Inf,
                     parallel.sz=1L,
                     mx.diff=TRUE,
                     tau=0.25,
                     ssgsea.norm=TRUE,
                     verbose=TRUE,
                     BPPARAM=SerialParam(progressbar=verbose)){

  kernelMethods <- match.arg(kernelMethods)

  ## filter genes according to verious criteria,
  ## e.g., constant expression
  exp <- .filterFeaturesgsea(exp)

  ## map to the actual features for which expression data is available
  mapped.pathwaylist <- .mapGeneSetsToFeatures(pathwaylist, rownames(exp))

  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user
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

  rval <- .ssgsea(exp, mapped.pathwaylist, kernelMethods, rnaseq, abs.ranking,
                  parallel.sz, mx.diff, tau, kernel, ssgsea.norm, verbose, BPPARAM)

  return(rval)

}



#' @keywords internal
#' @importFrom BiocParallel SerialParam
#' @importMethodsFrom BiocParallel bpworkers
#' @importFrom BiocParallel bpworkers
.ssgsea <- function(exp, pathwaylist,
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
    cat("Estimating UsessGSEA scores for", length(pathwaylist),"gene sets.\n")

  return(ssgsea(exp, pathwaylist, alpha=tau, parallel.sz=parallel.sz,
                normalization=ssgsea.norm, verbose=verbose, BPPARAM=BPPARAM))


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



#' @keywords internal
#' @importFrom DelayedArray t
#' @importFrom stats ecdf
#' @useDynLib PathActMarker
compute.gene.density <- function(exp, sample.idxs, rnaseq=FALSE, kernel=TRUE){
  n.test.samples <- ncol(exp)
  n.genes <- nrow(exp)
  n.density.samples <- length(sample.idxs)

  gene.density <- NA
  if (kernel) {
    # A = .C("matrix_density_R",
    #        as.double(t(exp[ ,sample.idxs, drop=FALSE])),
    #        as.double(t(exp)),
    #        R = double(n.test.samples * n.genes),
    #        n.density.samples,
    #        n.test.samples,
    #        n.genes,
    #        as.integer(rnaseq))$R
    A = double(n.test.samples * n.genes)
    matrix_density_R(
           as.double(t(exp[ ,sample.idxs, drop=FALSE])),
           as.double(t(exp)),
           A,
           n.density.samples,
           n.test.samples,
           n.genes,
           as.integer(rnaseq))

    gene.density <- t(matrix(A, n.test.samples, n.genes))
  } else {
    gene.density <- t(apply(exp, 1, function(x, sample.idxs) {
      f <- ecdf(x[sample.idxs])
      f(x)
    }, sample.idxs))
    gene.density <- log(gene.density / (1-gene.density))
  }

  return(gene.density)
}



#' @keywords internal
#' @importMethodsFrom BiocParallel bpiterate
#' @importMethodsFrom BiocParallel bplapply
#' @importFrom parallel splitIndices
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocParallel multicoreWorkers
compute.geneset.es <- function(exp, pathwaylist, sample.idxs, rnaseq=FALSE,
                               abs.ranking, parallel.sz=1L,
                               mx.diff=TRUE, tau=1, kernel=TRUE,
                               verbose=TRUE, BPPARAM=SerialParam(progressbar=verbose)) {
  num_genes <- nrow(exp)
  if (verbose) {
    if (kernel) {
      if (rnaseq)
        cat("Estimating ECDFs with Poisson kernels\n")
      else
        cat("Estimating ECDFs with Gaussian kernels\n")
    } else
      cat("Estimating ECDFs directly\n")
  }

  if (parallel.sz > 1 && length(sample.idxs > 100) && nrow(exp) > 100) {
    if (verbose)
      cat(sprintf("Estimating ECDFs in parallel on %d cores\n", as.integer(parallel.sz)))
    iter <- function(Y, n_chunks=BiocParallel::multicoreWorkers()) {
      idx <- splitIndices(nrow(Y), min(nrow(Y), n_chunks))
      i <- 0L
      function() {
        if (i == length(idx))
          return(NULL)
        i <<- i + 1L
        Y[idx[[i]], , drop=FALSE]
      }
    }
    gene.density <- bpiterate(iter(exp, 100),
                              compute.gene.density,
                              sample.idxs=sample.idxs,
                              rnaseq=rnaseq, kernel=kernel,
                              REDUCE=rbind, reduce.in.order=TRUE,
                              BPPARAM=BPPARAM)
  } else
    gene.density <- compute.gene.density(exp, sample.idxs, rnaseq, kernel)

  compute_rank_score <- function(sort_idx_vec){
    tmp <- rep(0, num_genes)
    tmp[sort_idx_vec] <- abs(seq(from=num_genes,to=1) - num_genes/2)
    return (tmp)
  }

  rank.scores <- rep(0, num_genes)
  sort.sgn.idxs <- apply(gene.density, 2, order, decreasing=TRUE) # n.genes * n.samples

  rank.scores <- apply(sort.sgn.idxs, 2, compute_rank_score)

  m <- bplapply(pathwaylist, ks_test_m,
                gene.density=rank.scores,
                sort.idxs=sort.sgn.idxs,
                mx.diff=mx.diff, abs.ranking=abs.ranking,
                tau=tau, verbose=verbose,
                BPPARAM=BPPARAM)
  m <- do.call("rbind", m)
  colnames(m) <- colnames(exp)

  return (m)
}



#' @keywords internal
ks_test_m <- function(gset_idxs, gene.density, sort.idxs, mx.diff=TRUE,
                      abs.ranking=FALSE, tau=1, verbose=TRUE){

  n.genes <- nrow(gene.density)
  n.samples <- ncol(gene.density)
  n.geneset <- length(gset_idxs)

  geneset.sample.es = double(n.samples)
  ks_matrix_R(as.double(gene.density),
              geneset.sample.es,
              as.integer(sort.idxs),
              n.genes,
              as.integer(gset_idxs),
              n.geneset,
              as.double(tau),
              n.samples,
              as.integer(mx.diff),
              as.integer(abs.ranking))


  return(geneset.sample.es)
}



#' @keywords internal
.rndWalk <- function(gSetIdx, geneRanking, j, R, alpha) {
  indicatorFunInsideGeneSet <- match(geneRanking, gSetIdx)
  indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
  indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
  stepCDFinGeneSet <- cumsum((abs(R[geneRanking, j])^alpha *
                                indicatorFunInsideGeneSet)) /
    sum((abs(R[geneRanking, j])^alpha *
           indicatorFunInsideGeneSet))
  stepCDFoutGeneSet <- cumsum(!indicatorFunInsideGeneSet) /
    sum(!indicatorFunInsideGeneSet)
  walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

  sum(walkStat)
}




#' @keywords internal
.fastRndWalk <- function(gSetIdx, geneRanking, j, Ra) {
  n <- length(geneRanking)
  k <- length(gSetIdx)
  idxs <- sort.int(match(gSetIdx, geneRanking))

  stepCDFinGeneSet2 <-
    sum(Ra[geneRanking[idxs], j] * (n - idxs + 1)) /
    sum((Ra[geneRanking[idxs], j]))


  stepCDFoutGeneSet2 <- (n * (n + 1) / 2 - sum(n - idxs + 1)) / (n - k)

  walkStat <- stepCDFinGeneSet2 - stepCDFoutGeneSet2

  walkStat
}




#' @keywords internal
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom BiocParallel SerialParam
#' @importFrom sparseMatrixStats colRanks
#' @importFrom DelayedArray t
ssgsea <- function(X, geneSets, alpha=0.25, parallel.sz,
                   normalization=TRUE, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {

  n <- ncol(X)

  print("Calculating ranks...")

  R <- t(sparseMatrixStats::colRanks(X, ties.method = "average"))
  mode(R) <- "integer"

  print("Calculating absolute values from ranks...")

  Ra <- abs(R)^alpha

  es <- bplapply(as.list(1:n), function(j) {
    geneRanking <- order(R[, j], decreasing=TRUE)
    es_sample <- lapply(geneSets, .fastRndWalk, geneRanking, j, Ra)

    unlist(es_sample)
  }, BPPARAM=BPPARAM)

  es <- do.call("cbind", es)


  if (normalization) {
    print("Normalizing...")

    es <- es[, 1:n, drop=FALSE] / (range(es)[2] - range(es)[1])
  }

  if (length(geneSets) == 1)
    es <- matrix(es, nrow=1)

  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)

  if(is(X, "dgCMatrix")){
    es <- as(as(as(es, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }

  return(es)
}






#' @keywords internal
#' @import methods
setGeneric("filterGeneSets", function(gSets, ...) standardGeneric("filterGeneSets"))



#' @keywords internal
#' @import methods
setMethod("filterGeneSets", signature(gSets="list"),
          function(gSets, min.sz=1, max.sz=Inf) {
            gSetsLen <- lengths(gSets)
            return (gSets[gSetsLen >= min.sz & gSetsLen <= max.sz])
          })



#' @keywords internal
#' @import methods
#' @importClassesFrom GSEABase GeneSetCollection
#' @importMethodsFrom GSEABase geneIds
setMethod("filterGeneSets", signature(gSets="GeneSetCollection"),
          function(gSets, min.sz=1, max.sz=Inf) {
            filterGeneSets(geneIds(gSets), min.sz, max.sz)
          })

















