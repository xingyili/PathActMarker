# Pathway activation evaluation algorithm：PLIER




#' returns the row names in common
#' @param exp1 One matrix with gene rownames
#' @param exp2 Another matrix with gene rownames
#' @keywords internal
commonRows=function(exp1, exp2){
  intersect(rownames(exp1), rownames(exp2))
}


#' z-score each row of the exp
#' @param x gene expression matrix, with genes in rows
#' @keywords internal
rowNorm=function(x){
  s=apply(x,1,sd)
  m=apply(x,1,mean);
  x=sweep(x,1,m)
  x=sweep(x,1,s,"/")
  return(x)
}


#' @keywords internal
pinv_ridge=function (m, alpha = 0){
  msvd = svd(m)
  if (length(msvd$d) == 0) {
    return(array(0, dim(m)[2:1]))
  }
  else {
    if (alpha > 0) {
      ss = (msvd$d^2) + alpha^2
      msvd$d = ss/msvd$d
    }
    out = msvd$v %*% (1/msvd$d * t(msvd$u))
    rownames(out) = rownames(m)
    colnames(out) = colnames(m)
    return(out)
  }
}

#' estimates the number of 'significant' principle components for the SVD decomposition: the minimum k for PLIER
#' @param  exp the same exp as to be used for PLIER (z-score recommended)
#' @param method Either "eblow" (fast) or "permutation" (slower, but less heuristic)
#' @param B number of permutations
#' @param seed seed for reproducibility
#' @import rsvd
#' @keywords internal
num_pc = function (exp,  method = "elbow", B = 20,  seed = NULL) {

  method = match.arg(method, c("elbow", "permutation"))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  warn <- NULL

  if ((class(exp) != "list") & (class(exp) != "rsvd")) {
    message("Computing svd")
    n <- ncol(exp)
    m <- nrow(exp)
    exp = rowNorm(exp)
    if (n < 500) {
      k = n
    }
    else{
      k = max(200, n / 4)
    }
    if (k == n) {
      temp <- svd(exp)
    }
    else{
      set.seed(123456)
      temp <- rsvd::rsvd(exp, k, q = 3)
    }
  }
  else if (!is.null(exp[["d"]])) {
    if (method == "permutation") {
      message("Original exp is needed for permutation method.\nSetting method to elbow")
      method = "elbow"
    }
    temp = exp
  }

  if (method == "permutation") {
    #nn = min(c(n, m))
    dstat <- temp$d[1:k] ^ 2 / sum(temp$d[1:k] ^ 2)
    dstat0 <- matrix(0, nrow = B, ncol = k)
    for (i in 1:B) {
      dat0 <- t(apply(exp, 1, sample, replace = FALSE))
      if (k == n) {
        uu0 <- svd(dat0)
      }
      else{
        set.seed(123456)

        uu0 <- rsvd:rsvd(dat0, k, q = 3)
      }
      dstat0[i,] <- uu0$d[1:k] ^ 2 / sum(uu0$d[1:k] ^ 2)
    }
    psv <- rep(1, k)
    for (i in 1:k) {
      psv[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:k) {
      psv[i] <- max(psv[(i - 1)], psv[i])
    }

    nsv <- sum(psv <= 0.1)
  }
  else if (method == "elbow") {
    x = smooth(xraw <- abs(diff(diff(temp$d))), twiceit = T)
    nsv = which(x <= quantile(x, 0.5))[1] + 1

  }
  return(nsv)
}


#'Solves for the U coefficients making efficient utilizatoin of the lasso path
#' @keywords  internal
#' @param Z current Z estimate
#' @param Chat the inverse of the C matrix
#' @param pathwaylist the prior pathway
#' @param penalty_factor Penalties for different pathways, must have size ncol(pathwaylist).
#' @param pathwaySelection Method to use for pathway selection.
#' @param glm_alpha The elsatic net alpha parameter
#' @param maxPath The maximum number of pathways to consider
#' @param target_frac The target fraction on non-zero columns
#' @param L3 Solve with a given L3, no search
#' @import glmnet
solveU=function(Z,  Chat, pathwaylist, penalty_factor,pathwaySelection="fast", glm_alpha=0.9, maxPath=10, target_frac=0.7, L3=NULL){


  Ur=Chat%*%Z
  Ur=apply(-Ur,2,rank)
  Urm=apply(Ur,1,min)

  U=matrix(0,nrow=ncol(pathwaylist), ncol=ncol(Z))
  if(is.null(L3)){

    lambdas=exp(seq(-4,-12,-0.125))
    results=list()
    lMat=matrix(nrow=length(lambdas), ncol=ncol(Z))
    for(i in 1:ncol(Z)){
      if(pathwaySelection=="fast"){
        iip=which(Ur[,i]<=maxPath)
      }else{
        iip=which(Urm<=maxPath)
      }#else
      gres=glmnet::glmnet(y=Z[,i], x=pathwaylist[,iip], penalty_factor = penalty_factor[iip], alpha=glm_alpha, lower.limits=0, lambda = lambdas,intercept=T,  standardize=F )

      gres$iip=iip
      lMat[,i]=colSums(as.matrix(gres$beta)>0)
      results[[i]]=gres
    }
    fracs=rowMeans(lMat>0)
    iibest=which.min(abs(target_frac-fracs))
    iibest


    for(i in 1:ncol(Z)){
      U[results[[i]]$iip,i]=results[[i]]$beta[,iibest]
    }#for i
    rownames(U)=colnames(pathwaylist)
    colnames(U)=1:ncol(Z)

    Utmp=solveU(Z,  Chat, pathwaylist, penalty_factor,pathwaySelection="fast", glm_alpha=0.9, maxPath=10,  L3=lambdas[iibest])

    #stop()
    return(list(U=U, L3=lambdas[iibest]))
  }
  else{
    for(i in 1:ncol(Z)){
      if(pathwaySelection=="fast"){
        iip=which(Ur[,i]<=maxPath)
      }else{
        iip=which(Urm<=maxPath)
      }#else
      gres=glmnet::glmnet(y=Z[,i], x=pathwaylist[,iip], penalty_factor = penalty_factor[iip], alpha=glm_alpha, lower.limits=0, lambda = L3,intercept=T,  standardize=F )
      U[iip,i]=as.numeric(gres$beta)
    }

    return(U)
  }
}



#' @keywords  internal
wrapString=function(string, width=30){

  string=lapply(string, function(s){ paste(strwrap(gsub("_", " ",s), width=width), collapse="\n")})
  unlist(string)
}


#' @keywords  internal
QV=function(pval){

  x=try(qvalue(pval))

  if(!is.list(x)){
    warning("Q-value error, defaulting to BH")
    #hist(pval)
    return(p.adjust(pval, method="BH"))
  }
  else{
    return(x$qvalue)
  }
}


#' @keywords  internal
mydist = function(x) {
  as.dist(1 - t(cor(t(x))))
}


#' @keywords  internal
BH= function(pval){p.adjust(pval, method="BH")}

#' SVD based smoothing for single cell RNAseq data
#' @param svdres svd result
#' @param k number of components to use
#' @keywords internal
DataSmooth=function(svdres,k){
  k=1:k
  ds=sweep(svdres$u[, k],2,svdres$d[k],"*")%*%t(svdres$v[, k])
  return(ds)
}

#' Rename pathway matrix gene names.
#' @param pathway the pathway.
#' @param map Gene name map.
#' @keywords internal
mapPathway=function(pathway, map){
  cm=commonRows(map, pathway)
  show(length(cm))
  pathway=pathway[cm,]
  rownames(pathway)=map[cm,1]
  return(pathway)
}


#' Creates a binary cell-type marker matrix using prior results.
#' @param plierResults A PLIER result
#' @param pathwaylist the binary prior information matrix that was used to compute the plierResult.
#' @param num The number of marker genes.
#' @param index The indecies of PLIER latent variables that are believed to represent cell-type proportions.
#' @keywords internal
plierResToMarkers=function(plierResults, pathwaylist, num=20, index=NULL){

  ii=which(colSums(plierResults$U)>0)
  if(! is.null(index)){
    ii=intersect(ii,index)
  }
  Zuse=plierResults$Z[,ii, drop=F]

  for(i in 1:length(ii)){
    lv=ii[i]
    paths=names(which(plierResults$U[,lv]<0.01))
    genes=names(which(rowSums(pathwaylist[,paths])>0))
    genesNotInPath=setdiff(rownames(Zuse), genes)
    Zuse[genesNotInPath,i]=0
  }

  tag_matrix=apply(-Zuse,2,rank)
  colnames(tag_matrix)=rownames(plierResults$B)[ii]
  iim=apply(tag_matrix,1,min)
  iig=which(iim<=num)
  tag_matrix=tag_matrix[iig,]
  iin=rowSums(tag_matrix<=num)
  iimulti=which(iin>1)
  if(length(iimulti)>0){
    message(paste0("Genes not matched uniquely: ", paste(names(iimulti), collapse=", ")))
  }
  tag_matrix=(tag_matrix<=num)+1-1

  return(tag_matrix)
}



#' @keywords internal
AUC<-function(labels, values){
  posii=which(labels>0)
  negii=which(labels<=0)
  posn=length(posii)
  negn=length(negii)
  posval=values[posii]
  negval=values[negii]
  myres=list()
  if(posn>0&negn>0){
    res=wilcox.test(posval, negval, alternative="greater", conf.int=TRUE);

    myres$low=res$conf.int[1]
    myres$high=res$conf.int[2]
    myres$auc=(res$statistic)/(posn*negn)
    myres$pval=res$p.value
  }
  else{
    myres$auc=0.5
    myres$pval=NA
  }
  return(myres)
}


#' get the p-value cutoff for a specific FDR
#'
#' @param plierResults A PLIER result
#' @param fdr.cutoff The cross-validation significance cutoff for a pathway.
#' @keywords internal
getCutoff=function(plierResults,  fdr.cutoff=0.01){
  max(plierResults$summary[plierResults$summary[,"FDR"]<=fdr.cutoff,"p-value"])
}


#' names latent variables according to their pathway useage
#'
#' @param plierResults A PLIER result
#' @param top The number of pathway to use.
#' @param fdr.cutoff The cross-validation significance cutoff for a pathway to be considered for naming. If no pathways satisfy the cutoff the raw coefficients are used.
#' @param use of coef or AUC, whether LVs are named based on U coefficients or AUCs. Defualt: coef.
#' @keywords internal
nameB=function(plierResults, top=1, fdr.cutoff=0.01, use=c("coef", "AUC")){
  use=match.arg(use, c("coef", "AUC"))
  names=vector("character",ncol(plierResults$U))
  if(use=="coef"){
    Uuse=plierResults$U
  }
  else{
    Uuse=plierResults$Uauc
  }
  if(!is.null(plierResults[["Up"]])){
    pval.cutoff=max(plierResults$summary[plierResults$summary[,5]<fdr.cutoff,4])

    Uuse[plierResults$Up>pval.cutoff]=0

  }
  else{
    warning("No p-values in PLIER object: using coefficients only")
  }
  mm=apply(Uuse,2,max)
  for(i in 1:ncol(plierResults$U)){
    if(mm[i]>0){
      names[i]=paste(i,names(sort(Uuse[,i],T))[1:top], sep=",")
    }
    else if(max(plierResults$U[,i])>0){
      names[i]=paste(i,names(sort(plierResults$U[,i],T))[1:top], sep=",")
    }
    else{
      names[i]=paste("LV",i)
    }
  }

  return(names)

}

#' Computes the ridge pseudo-inverse of the prior information matrix.
#' @param gsMat The prior information matrix. The genes have to match the gene expression exp.
#' @param lambda The regularization paramter
#' @keywords internal
computeChat=function(gsMat, lambda=5){
  Chat=pinv_ridge(crossprod(gsMat,), lambda)%*%(t(gsMat))
}


#' @keywords internal
copyMat=function(mat, zero=F){
  matnew=matrix(nrow=nrow(mat), ncol=ncol(mat))
  rownames(matnew)=rownames(mat)
  colnames(matnew)=colnames(mat)
  if(zero)
    matnew[]=0
  return(matnew)
}


#' crossVal
#' @keywords internal
#' @param pathwaylist the real prior info matrix
#' @param pathwaylistcv the zeroed-out prior info matrix used for PLIER computations
crossVal=function(plierResults, exp, pathwaylist, pathwaylistcv){

  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(plierResults$U)>0)
  Uauc=copyMat(plierResults$U,T)
  Up=copyMat(plierResults$U,T)
  Up[]=1
  for ( i in ii){

    iipath=which(plierResults$U[,i]>0)

    if (length(iipath) > 1){
      for(j in iipath){
        iiheldout=which((rowSums(pathwaylist[,iipath, drop=F])==0) |(pathwaylist[,j]>0&pathwaylistcv[,j]==0))
        aucres=AUC(pathwaylist[iiheldout,j], plierResults$Z[iiheldout,i])
        out=rbind(out,c(colnames(pathwaylist)[j], i, aucres$auc, aucres$pval))
        Uauc[j,i]=aucres$auc
        Up[j,i]=aucres$pval
      }}else{
        j <- iipath
        iiheldout=which((rowSums(matrix(pathwaylist[,iipath],ncol=1))==0) |(pathwaylist[,j]>0&pathwaylistcv[,j]==0))
        aucres=AUC(pathwaylist[iiheldout,j], plierResults$Z[iiheldout,i])
        out=rbind(out,c(colnames(pathwaylist)[j], i, aucres$auc, aucres$pval))
        Uauc[j,i]=aucres$auc
        Up[j,i]=aucres$pval
      }#else
  }
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5]=BH(out[,4])
  colnames(out)=c("pathway", "LV index", "AUC", "p-value", "FDR")
  return(list(Uauc=Uauc, Upval=Up, summary=out))
}

#' getAUC
#' @keywords internal
#' @param plierResults current PLIER result
#' @param exp the input exp
#' @param pathwaylist the  prior info matrix
getAUC=function(plierResults, exp, pathwaylist){
  Y=exp
  B=plierResults$B
  Z=plierResults$Z
  Zcv=copyMat(Z)
  k=ncol(Z)
  L1=plierResults$L1
  L2=plierResults$L2
  for (i in 1:5){
    ii=(0:(floor(nrow(exp)/5)-1))*5+i
    ii=ii[ii<=nrow(Z)]


    Bcv=solve(crossprod(Z[-ii,])+L2*diag(k))%*%t(Z[-ii,])%*%Y[-ii,]

    Zcv[ii,]=Y[ii, ]%*%t(Bcv)%*%solve(tcrossprod(Bcv)+L1*diag(k))
  }

  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(plierResults$U)>0)
  Uauc=copyMat(plierResults$U,T)
  Up=copyMat(plierResults$U,T)
  Up[]=1;
  for ( i in ii){

    iipath=which(plierResults$U[,i]>0)

    for(j in iipath){
      aucres=AUC(pathwaylist[,j], Zcv[,i])
      out=rbind(out,c(colnames(pathwaylist)[j], i, aucres$auc, aucres$pval))
      Uauc[j,i]=aucres$auc
      Up[j,i]=aucres$pval
    }
  }
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5]=BH(out[,4])
  colnames(out)=c("pathway", "LV index", "AUC", "p-value", "FDR")

  return(list(Uauc=Uauc, Upval=Up, summary=out))
}


#' @keywords internal
nonEstimable=function (x)
{
  x = as.matrix(x)
  p = ncol(x)
  QR = qr(x)
  if (QR$rank < p) {
    n = colnames(x)
    if (is.null(n))
      n = as.character(1:p)
    notest = n[QR$pivot[(QR$rank + 1):p]]
    blank = notest == ""
    if (any(blank))
      notest[blank] = as.character(((QR$rank + 1):p)[blank])
    return(notest)
  }
  else {
    return(NULL)
  }
}


#' @keywords internal
resid=function(dat, lab, useMean=T){
  if (is.null(dim(lab))){
    mod=model.matrix(~1+lab);
  }
  else{
    mod=lab
  }
  ne = nonEstimable(mod)
  if (!is.null(ne)){
    cat("Coefficients not estimable:", paste(ne, collapse = " "),
        "\n")
    mod=mod[, -match(ne, colnames(mod))]
  }

  n=dim(dat)[2]
  Id=diag(n)
  out=dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*%
                 t(mod))
  colnames(out)=colnames(dat)
  if (useMean){
    out=sweep(out,1,apply(dat,1,mean), "+")
  }

  return(out)
}




#' @keywords internal
scad=function(x, lambda,a=3.7){

  iip=which(abs(x)>2*lambda & abs(x)<a*lambda)

  iin=which(abs(x)<=2*lambda)

  x[iin]=pmax(0, abs(x[iin])-lambda)*sign(x[iin])
  x[iip]=((a-1)*x[iip]-sign(x[iip])*a*lambda)/(a-2)

  return(x)
}


#' @keywords internal
quicksoft=function (x, d) {
  return(sign(x) * pmax(0, abs(x) - d))
}


#' @keywords internal
scadZ=function(Z, ii=1:ncol(Z), lambda){
  Zn=colSumNorm(Z, return.all = T)
  Zt=Z

  Zt[,ii]=apply(Z[,ii], 2, function(x){scad(x, lambda)})
  return(Zt)
}


#' @keywords internal
softZ=function(Z, ii=1:ncol(Z), lambda){
  Zn=colSumNorm(Z, return.all = T)
  Zt=Z

  Zt[,ii]=apply(Z[,ii], 2, function(x){quicksoft(x, lambda)})
  return(Zt)
}


#' @keywords internal
getEnrichmentVals=function(plierResults, pathwayMat, ngenes=50,auc.cutoff=0.7, fdr.cutoff=0.01){
  pathwayMat=pathwayMat[rownames(plierResults$Z), rownames(plierResults$U)]
  Uuse=plierResults$U
  Uuse[plierResults$Uauc<auc.cutoff]=0
  Uuse[plierResults$Up>getCutoff(plierResults, fdr.cutoff)]=0
  inpath=intop=double(ncol(plierResults$Z))

  for(i in 1:ncol(plierResults$Z)){
    iipath=which(Uuse[,i]>0)
    if(length(iipath)>0){
      pathGenes=names(which(rowSums(pathwayMat[,iipath, drop=F])>0))
      topGenes=names(sort(plierResults$Z[,i], T)[1:ngenes])
      pathGenesInt=intersect(pathGenes, topGenes)
      inpath[i]=length(pathGenes)
      intop[i]=length(pathGenesInt)
    }}
  return(cbind(intop/inpath,intop, inpath))
}


#' @keywords internal
tscale=function(x, zeroNA=T){
  s = apply(x, 1, sd, na.rm=T)
  m = apply(x, 1, mean, na.rm=T)
  x = sweep(x, 1, m)
  x = sweep(x, 1, s, "/")
  if(zeroNA){
    x[is.na(x)]=0
  }
  return(x)
}



#' Main PLIER function
#'
#' @param exp the exp to be processed with genes in rows and samples in columns. Should be z-scored or set scale=T
#' @param pathwaylist the  prior information
#' @param summary Display LV and p-value corresponding to different pathway
#' @param svdres Pre-computed result of the svd decomposition for exp
#' @param k The number of latent variables to return, leave as NULL to be set automatically using the num_pc "elbow" method
#' @param L1 L1 constant, leave as NULL to automatically select a value
#' @param L2 L2 constant, leave as NULL to automatically select a value
#' @param L3 L3 constant, leave as NULL to automatically select a value. Sparsity in U should be instead controlled by setting frac
#' @param frac The fraction of LVs that should have at least 1 prior inforamtion association, used to automatically set L3
#' @param max.iter Maximum number of iterations to perform
#' @param trace Display progress information
#' @param scale Z-score the exp before processing
#' @param Chat A ridge inverse of pathwaylist, used to select active pathways, expensive to compute so can be precomputed when running PLIER multiple times
#' @param maxPath The maximum number of active pathways per latent variable
#' @param doCrossval Whether or not to do real cross-validation with held-out pathway genes. Alternatively, all gene annotations are used and only pseudo-crossvalidation is done. The latter option may be preferable if some pathways of interest have few genes.
#' @param penalty_factor A vector equal to the number of columns in pathwaylist. Sets relative penalties for different pathway/geneset subsets. Lower penalties will make a pathway more likely to be used. Only the relative values matter. Internally rescaled.
#' @param glm_alpha Set the alpha for elastic-net
#' @param minGenes The minimum number of genes a pathway must have to be considered
#' @param tol Convergence threshold
#' @param seed Set the seed for pathway cross-validation
#' @param  allGenes Use all genes. By default only genes in the pathwaylist matrix are used.
#' @param rseed Set this option to use a random initialization, instead of SVD
#' @param pathwaySelection Pathways to be optimized with elstic-net penalty are preselected based on ridge regression results. "Complete" uses all top  pathways to fit individual LVs. "Fast" uses only the top pathways for the single LV in question.
#' @export
#' @examples
#' data(KICH)
#' data(KEGGgenesetsID)
#' KICHfilter <- NormalizeExp(KICH)
#' PLIERResults = UsePLIER(KICHfilter, KEGGgenesetsID)
UsePLIER=function(exp, pathwaylist,summary=F,svdres=NULL, k=NULL, L1=NULL, L2=NULL, L3=NULL,  frac=0.7,  max.iter=350, trace=F, scale=T, Chat=NULL, maxPath=10, doCrossval=T, penalty_factor=rep(1,ncol(pathwaylist)), glm_alpha=0.9, minGenes=10, tol=1e-6, seed=123456, allGenes=F, rseed=NULL, pathwaySelection=c("complete", "fast")){

  pathwaySelection=match.arg(pathwaySelection, c("complete", "fast"))

  genes <- unique(unlist(pathwaylist))
  adj_matrix <- matrix(0, nrow = length(genes), ncol = length(pathwaylist),
                       dimnames = list(genes, names(pathwaylist)))

  for (i in seq_along(pathwaylist)) {
    for (gene in pathwaylist[[i]]) {
      adj_matrix[gene, names(pathwaylist)[i]] <- 1
    }
  }
  pathwaylist <- adj_matrix

  if(scale){
    Y=rowNorm(exp)
  }else{
    Y=exp
  }

  if(nrow(pathwaylist)!=nrow(exp) || !all(rownames(pathwaylist)==rownames(exp))){
    if(!allGenes){
      cm=commonRows(exp, pathwaylist)
      message(paste("Selecting common genes:", length(cm)))
      pathwaylist=pathwaylist[cm,]
      Y=Y[cm,]
    }
    else{
      extra.genes=setdiff(rownames(exp), rownames(pathwaylist))
      eMat=matrix(0, nrow=length(extra.genes), ncol=ncol(pathwaylist))
      rownames(eMat)=extra.genes
      pathwaylist=rbind(pathwaylist, eMat)
      pathwaylist=pathwaylist[rownames(exp),]
    }

  }

  numGenes=colSums(pathwaylist)
  heldOutGenes=list()
  iibad=which(numGenes<minGenes)
  pathwaylist[, iibad]=0
  message(paste("Removing", length(iibad), "pathways with too few genes"))

  if(doCrossval){

    pathwaylistCV=pathwaylist
    if(!is.null(seed))
      set.seed(seed)
    for(j in 1:ncol(pathwaylistCV)){

      iipos=which(pathwaylistCV[,j]>0)
      iiposs=sample(iipos, length(iipos)/5)
      pathwaylistCV[iiposs,j]=0
      heldOutGenes[[colnames(pathwaylist)[j]]]=rownames(pathwaylist)[iiposs]
    }
    C = pathwaylistCV
  }else{
    C=pathwaylist
  }


  nc=ncol(pathwaylist)
  ng=nrow(exp)
  ns=ncol(exp)

  Bdiff=-1
  BdiffTrace=double()
  BdiffCount=0
  if(is.null(Chat)){
    Cp=crossprod(C)
    Chat=pinv_ridge(crossprod(C), 5)%*%(t(C))
  }
  YsqSum=sum(Y^2)

  if(!is.null(svdres) && nrow(svdres$v)!=ncol(Y)){
    message("SVD V has the wrong number of columns")
    svdres=NULL
  }

  if(is.null(svdres)){
    message("Computing SVD")
    if(ns>500){
      message("Using rsvd")
      set.seed(123456);svdres=rsvd::rsvd(Y, k=min(ns, max(200, ns/4)), q=3)
    }
    else{
      svdres=svd(Y)
    }
    message("Done")
  }

  if(is.null(k)){
    k=num_pc(svdres)*2
    k <- min(k, floor(ncol(Y)*0.9))
    message("k is set to ", k)
  }


  if(is.null(L2)){
    show(svdres$d[k])
    L2=svdres$d[k]
    print(paste0("L2 is set to ",L2))
  }
  if(is.null(L1)){
    L1=L2/2
    print(paste0("L1 is set to ",L1))
  }


  B=t(svdres$v[1:ncol(Y), 1:k]%*%diag(svdres$d[1:k]))
  Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
  Z[Z<0]=0
  if(!is.null(rseed)){
    message("using random start")
    set.seed(rseed)
    B=t(apply(B, 1, sample))
    Z=apply(Z,2,sample)
  }

  U=matrix(0,nrow=ncol(C), ncol=k)

  round2=function(x){signif(x,4)}
  message(paste0("errorY (SVD based:best possible) = ", round2(mean((Y-Z%*%B)^2))))

  iter.full.start=iter.full=20

  curfrac=0
  nposlast=Inf
  npos=-Inf
  if(!is.null(L3)){
    L3.given=T
  }else{
    L3.given=F
  }

  for ( i in 1:max.iter){

    if(i>=iter.full.start){

      if(i==iter.full & !L3.given){
        Ulist=solveU(Z, Chat, C, penalty_factor, pathwaySelection, glm_alpha, maxPath, target_frac = frac)
        U=Ulist$U
        L3=Ulist$L3
        message(paste("New L3 is", L3))
        iter.full=iter.full+iter.full.start
      }
      else{

        U=solveU(Z, Chat, C, penalty_factor, pathwaySelection, glm_alpha, maxPath, L3=L3)
      }
      curfrac=(npos<-sum(apply(U,2,max)>0))/k
      Z1=Y%*%t(B)
      Z2=L1*C%*%U
      ratio=median((Z2/Z1)[Z2>0&Z1>0])
      Z=(Z1+Z2)%*%solve(tcrossprod(B)+L1*diag(k))
    }

    else{
      Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
    }

    Z[Z<0]=0
    oldB=B
    B=solve(t(Z)%*%Z+L2*diag(k))%*%t(Z)%*%Y

    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)

    err0=sum((Y-Z%*%B)^2)+sum((Z-C%*%U)^2)*L1+sum(B^2)*L2
    if(trace & i >=iter.full.start){

      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", prior information ratio= ", round(ratio,2), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))), ";pos. col. U=", sum(colSums(U)>0))
    }
    else if (trace){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }

    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
      message("Bdiff is not decreasing")
    }
    else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }

    if(Bdiff<tol){
      message(paste0("converged at  iteration ", i))
      break
    }
    if( BdiffCount>5){
      message(paste0("converged at  iteration ", i, " Bdiff is not decreasing"))
      break
    }

  }

  rownames(U)=colnames(pathwaylist)
  colnames(U)=rownames(B)=paste0("LV", 1:k)
  out=list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C, L1=L1, L2=L2, L3=L3, heldOutGenes=heldOutGenes)

  if(doCrossval){
    outAUC=crossVal(out, Y, pathwaylist, pathwaylistCV)
  }
  else{
    message("Not using cross-validation. AUCs and p-values may be over-optimistic")
    outAUC=getAUC(out, Y, pathwaylist)
  }
  out$withPrior=which(colSums(out$U)>0)
  out$Uauc=outAUC$Uauc
  out$Up=outAUC$Upval
  out$summary=outAUC$summary
  tt=apply(out$Uauc,2,max)
  message(paste("There are", sum(tt>0.70), " LVs with AUC>0.70"))

  rownames(out$B)=nameB(out)

  if(summary){
    return(out$summary)
  }else{
    return(out$B)
  }
}

