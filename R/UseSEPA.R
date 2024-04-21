# Pathway activation evaluation algorithm：SEPA


#' @keywords internal
getData <- function(exp, SparseAdj){
  geneList <- rownames(exp)
  delIdx <- c()
  nodes1 <- SparseAdj[[1]]
  nodes2 <- SparseAdj[[2]]
  nodeSet <- c()
  for (i in 1:length(nodes1)){
    if (!(nodes1[i] %in% geneList)||
        !(nodes2[i] %in% geneList)){
      delIdx <- c(delIdx, i)
      next
    }
    if (!(nodes1[i] %in% nodeSet)){
      nodeSet <- c(nodeSet, nodes1[i])
    }
    if (!(nodes2[i] %in% nodeSet)){
      nodeSet <- c(nodeSet, nodes2[i])
    }
  }
  if (length(delIdx) != 0){
    nodes1 <- nodes1[-delIdx]
    nodes2 <- nodes2[-delIdx]
  }
  nodeNum = length(nodeSet)
  edgeNum = length(nodes1)
  singleAdj = matrix(0, nrow = nodeNum, ncol = nodeNum)
  rownames(singleAdj) <- nodeSet
  colnames(singleAdj) <- nodeSet
  for (i in 1:edgeNum){
    singleAdj[as.character(nodes1[i]), as.character(nodes2[i])] = 1
    singleAdj[as.character(nodes2[i]), as.character(nodes1[i])] = 1
  }
  sampleNum = ncol(exp)
  singleExp = c()
  for (i in nodeSet){
    singleExp = c(singleExp, exp[i,])
  }
  singleExp <- matrix(singleExp, ncol = sampleNum, byrow = TRUE)
  singleExp <- as.data.frame(singleExp)
  singleAdj <- as.data.frame(singleAdj)
  return(list(singleExp, singleAdj))
}

#' @keywords internal
matfind <- function (x)
{
  expr <- if (is.logical(x)) {
    x
  }
  else {
    x != 0
  }
  indices <- which(expr)
  col_indices <- (indices - 1) %/% nrow(expr) + 1
  row_indices <- (indices - 1) %% nrow(expr) + 1
  return(list(row_indices, col_indices))

}


#' Title Pathway activation evaluation algorithm：SEPA
#'
#' @param exp the gene expression matrix.
#' @param pathwaySparseAdj the pathway genesets.
#' @return a matrix with rows corresponding to the pathways and columns corresponding to the patients.
#' @export
#' @examples
#' data(KICH)
#' data(KEGGgenenetworkID)
#' KICHfilter <- NormalizeExp(KICH)
#' SEPAResults <- UseSEPA(KICHfilter, KEGGgenenetworkID)
# SEPA:main function
UseSEPA <- function(exp, pathwaySparseAdj){
  sampleNum = ncol(exp)
  ret = c()
  adjusted_exp <- exp
  for (i in 1:ncol(exp)) {
    adjusted_exp[,i] <- (exp[,i] - min(exp[,i])) / (max(exp[,i]) - min(exp[,i]))
  }
  for (adj in pathwaySparseAdj){
    dataInfo <- getData(adjusted_exp, adj)
    geneExpressionData <- dataInfo[[1]]
    adjMatrix <- dataInfo[[2]]
    esp <- 1e-6
    geneExpressionData <- geneExpressionData + esp

    geneNum <- nrow(geneExpressionData)
    geneSample <- length(geneExpressionData)

    res <- rep(0, geneSample)
    diag(adjMatrix) <- 1
    indices <- matfind(adjMatrix)
    row <- indices[[1]]
    col <- indices[[2]]

    stopError <- 10^-1
    maxIter <- 10^4
  
    for (i in 1:geneSample){
      singleSample <- geneExpressionData[, i]/sum(geneExpressionData[, i])
      varB <- Matrix::sparseMatrix(row, col, x = singleSample[row])
      varB <- as.matrix(varB)

      flagIter <- 1
      curIter <- 0
      alphaPre <- as.matrix(singleSample)
      betaPre <- matrix(1, nrow = geneNum)
      varA <- as.matrix(adjMatrix)

      while(flagIter){
        alpha <- 1./(varA%*%betaPre)
        beta <- singleSample / (t(varB) %*% alpha)
        curError <- max(abs(rbind(alpha-alphaPre, beta-betaPre)))
        alphaPre <- alpha
        betaPre <- beta
        curIter <- curIter + 1
        flagIter <- (curError > stopError && curIter < maxIter)
      }

      res[i] <- -((singleSample) %*% c(log(alpha * beta)) )
    }

    eigval <- (eigen(adjMatrix))$values
    modEigval <- Mod(eigval)
    lamda <- max(modEigval)

    maxSEPA <- log(lamda)
    res <- res/maxSEPA
    ret <- c(ret, res)
  }
  ret <- matrix(ret, ncol = sampleNum, byrow = TRUE)
  rownames(ret) <- names(pathwaySparseAdj)
  colnames(ret) <- colnames(exp)
  return (ret)
}



