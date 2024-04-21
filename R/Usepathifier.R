# Pathway activation evaluation algorithm：pathifier


#' @keywords internal
#' @importFrom  R.oo trim
.getmeasuredgenesinpathway = function(paths,allgenes) {

  len <-length(paths)
  pathways<-vector('list',len)
  for (i in 1:len) {
    n<-length(paths[[i]])
    isin<-matrix(FALSE,n)
    for(j in 1:n) {
      isin[j]<-(length(grep(paste('\\b',R.oo::trim(paths[[i]][j]),'\\b',sep=""),allgenes))>0)
    }
    pathways[[i]]<-unique(paths[[i]][isin])
  }
  return(pathways)
}


#' @keywords internal
#' @importFrom  R.oo trim
.getpathway = function(paths,allgenes,exp) {

  len<-length(paths)
  x<-NULL
  isin<-rep(FALSE,len);
  for(i in 1:len) {
    ind<-unique(grep(paste('\\b',R.oo::trim(paths[i]),'\\b',sep=""),allgenes))
    n<-length(ind)
    if (n>0) {
      if (n==1)
        t<-exp[ind,]
      else
        t<-colMeans(exp[ind,])
      if (var(t)>0) {
        isin[i]=TRUE;
        x<-c(x,t)
      }
    }
  }
  if (is.null(x)) {
    list(x=NULL,isin=isin)
  } else {
    list(x=matrix(x,nrow=ncol(exp)),isin=isin)
  }
}




#' @keywords internal
#' @importFrom  princurve principal_curve
.score_pathway =  function(x,m,ranks,calcerr=FALSE,thresh = 0.0005,maxit=200,start,logfile = "") {
  x<-x[,apply(x,2,sd)>0.001]
  k<-dim(x)[2]
  if (k<3) {
    c<-NULL
    cat(file=logfile,append=TRUE,'scoring failed (k=',k,').\n')
  } else {
    d<-matrix(0,1,m)
    if (start == "by pca") {
      start <- NULL
    } else if (start == "by ranks") {
      start <- aggregate(x, by=list(ranks), FUN=mean)
      start <- as.matrix(start[,-1])
    }
    c<-princurve::principal_curve(x, start=start, thresh=thresh, maxit=maxit, tol = 1e-10)
  }
  if (!is.null(c)) {
    d[c$ord[1]]=0
    for (j in 2:m) {
      d[c$ord[j]]<-d[c$ord[j-1]]+dist(c$s[c$ord[(j-1):j],])
    }
    d=d/d[c$ord[m]]
    if (calcerr) {
      e<-matrix(0,1,k)
      for (i in 1:k) {
        e[i]<-mean((c$s[,i]-x[,i])^2)
      }
    } else {
      e <- FALSE;
    }
    list(score=d,error=e,thecurve=c)
  } else {
    cat(file=logfile,append=TRUE,'scoring failed.\n')
    NULL
  }
}



#' @keywords internal
.samplings_stdev = function(m,n,attempts,z,ranks,samplings,start,logfile = "") {

  dall<-array(dim=c(attempts,n))
  skip<-0
  for(a in 1:attempts) {
    res<-.score_pathway(z[samplings[a,],],m,ranks[samplings[a,]],start=start,logfile=logfile)
    if (!is.null(res)) {
      dall[a,samplings[a,]] <- res$score
    } else {
      skip <- skip+1
    }
  }
  if (skip < attempts/2) {
    mean(apply(dall,2,sd,'na.rm'=TRUE), 'na.rm'=TRUE)
  } else {
    Inf
  }

}




#' @keywords internal
.score_all_pathways_helper = function(z, ranks, samplings, i, attempts, maximize_stability, logfile = "",start) {

  n<-dim(z)[1]
  k<-dim(z)[2]
  m<-dim(samplings)[2]
  mincheck<-5
  kmin=max(floor(0.8*k),mincheck+1)
  mindelta=min(0.009,max(0.002,1.5/k))
  sig<-matrix(0,1,k)
  res<-.score_pathway(z,n,ranks,calcerr=TRUE,start=start,logfile=logfile)
  if (is.null(res)) {
    cat(file=logfile,append=TRUE,'pathway ', i, '> scoring failed 1.\n')
  } else {
    sig<-.samplings_stdev(m,n,attempts,z,ranks,samplings,start=start)
    if (sig>10000) {
      cat(file=logfile,append=TRUE,'pathway ', i, '> scoring failed 2 (sig:', sig, ').\n')
      res<-NULL
    } else {
      origsig<-sig
      cat(file=logfile,append=TRUE,'pathway ', i, '> sig:', sig, '\n')
      isin<-1:k
      if (maximize_stability) {
        testsig<-max(mincheck,floor(0.1*k))
        newsig<-rep(0,testsig)
        while ((k>=kmin)&(sig>0.05)) {
          se<-sort(res$error,index.return=TRUE,decreasing=TRUE)
          for (j in 1:testsig) {
            newsig[j]<-.samplings_stdev(m,n,attempts,z[,-se$ix[j]],ranks,samplings,start=start)
          }
          wj<-which.min(newsig)
          cat(file=logfile,append=TRUE,'pathway ', i, ' k=', k, '(', ncol(res$thecurve$s), ') wj=', wj, '>new sig:', newsig[wj])
          if (sig-newsig[wj]<mindelta) {
            cat(file=logfile,append=TRUE,' x rejected\n')
            break
          }
          cat(file=logfile,append=TRUE,' | accepted!\n')
          sig<-newsig[wj]
          isin<-isin[-se$ix[wj]];
          z<-z[,-se$ix[wj]]
          k<-k-1
          res<-.score_pathway(z,n,ranks,calcerr=TRUE,start=start,logfile=logfile)
          if (is.null(res)) {
            cat(file=logfile,append=TRUE,'pathway ', i, '> scoring failed 3.\n')
            break;
          }
        }
      }
    }
  }
  if (is.null(res)) {
    NULL
  } else {
    list(score=res$score,thecurve=res$thecurve,z=z,isin=isin,sig=sig,origsig=origsig,k=k)
  }
}




#' Title pathifier,a pathway activation evaluation algorithm.
#'
#' @param exp gene expression matrix, where n is the number of genes and m the number of samples
#' @param paths a list of p pathways
#' @param label  normal(0) and tumor (1) in patient samples.
#' @param ranks external knowledge on the ranking of the m samples, if exists (to use initial guess)
#' @param attempts number of repeats  to stability
#' @param maximize_stability  throw away components leading to low stability of sampling noise.
#' @param logfile  name of the file the log should be written to (use stdout if empty).
#' @param samplings a matrix specifying the samples that should be chosen in each sampling attempt, chooses a random matrix if samplings is NULL.
#' @param min_exp the minimal expression considered as a real signal. Any values below are thresholded to be min_exp.
#' @param min_std the minimal allowed standard deviation of each gene.
#'
#' @return the iES_mat obtained by pathifier
#' @export
#'
#' @examples
#' data(KICH)
#' data(KEGGgenesetsID)
#' KICHfilter <- NormalizeExp(KICH)
#' pathifierResults = Usepathifier(KICHfilter, KEGGgenesetsID)
Usepathifier =  function(exp, paths, label = NULL, ranks = NULL, attempts = 100, maximize_stability = TRUE, logfile = "", samplings=NULL, min_exp=4, min_std=0.4) {
  allgenes<-rownames(exp)
  pathwaynames =  names(paths)
  # label <- as.logical(1 - as.numeric(label))
  if (!is.null(label))
  {
    label <- as.logical(1 - as.numeric(label))
  }
  cat(file=logfile,append=FALSE,'robust_score_bydist. min_exp=',min_exp,', min_std=',min_std,'\n')
  exp[exp<min_exp]=min_exp;
  n<-ncol(exp)
  if (is.null(label)) {
    label <- rep(TRUE,n);
    start <- "by pca";
  } else {
    start <- "by ranks"
  }
  if (is.null(ranks)) ranks <- !label;
  ranks <- rank(ranks)
  if ((length(label)!=n)||(length(ranks)!=n)) {
    stop("invalid dimentions");
  }
  len<-length(paths)
  nn<-sum(label)
  m<-floor(0.8*(n-nn))+nn
  if (is.null(samplings)) {
    samplings<-matrix(0,attempts,m)
    w<-which(!label)
    for(a in 1:attempts) {
      samplings[a,]<-sort(c(w[sample(n-nn,m-nn)],which(label)))
    }
  }
  s<-NULL
  ind<-NULL
  for (i in 1:len) {

    pathway<-paths[[i]]
    pathwayindata<-.getpathway(pathway,allgenes,exp)
    k1=sum(pathwayindata$isin)
    if (k1<3) {
      si<-NULL
      cat(file=logfile,append=TRUE,'skipping pathway ',i,' k1=', k1,'\n')
    } else {
      x<-pathwayindata$x
      pathway<-pathway[pathwayindata$isin]
      xm<-colMeans(x[label,])
      xs<-apply(x[label,],2,sd)
      xs[xs<min_std]=min_std;
      if (0 %in% xs) {
        si<-NULL
        cat(file=logfile,append=TRUE,'skipping pathway ',i,' (0 in xs)\n')
      } else {
        z<-(x-matrix(rep(xm,each=n),nrow=n))/(matrix(rep(xs,each=n),nrow=n))
        t<-prcomp(z)
        k2=max(sum(t$sdev>1.1),4)
        k2=min(k2,k1,0.75*dim(x)[1],sum(t$sdev>0.25))
        if (k2<3) {
          si<-NULL
          cat(file=logfile,append=TRUE,'skipping pathway ',i,' k2=', k2,'\n')
        } else {
          pca<-t$x[,1:k2]
          res<-.score_all_pathways_helper(pca, ranks, samplings, i, attempts, maximize_stability, logfile, start=start)
          if (is.null(res)) {
            si<-NULL
            cat(file=logfile,append=TRUE,'skipping pathway ',i,'\n')
          } else {
            ind<-c(ind,i)
            si<-list(res$score,pathway,res$sig,res$origsig,res$k,res$thecurve$s,res$thecurve$ord,res$z,res$isin,xm,xs,t$center,t$rotation,k2)
          }
        }
      }
    }
    s<-rbind(s,si)
  }
  cat(file=logfile,append=TRUE,length(ind),'pathways processed with start=',start,'\n')
  rownames(s)<-pathwaynames[ind]
  results = list(scores = s[,1], genesinpathway=s[,2], newmeanstd=s[,3], origmeanstd=s[,4], pathwaysize=s[,5], curves=s[,6], curves_order=s[,7], z=s[,8],compin=s[,9],xm=s[,10],xs=s[,11],center=s[,12],rot=s[,13],pctaken=s[,14],samplings=samplings,sucess=ind,logfile=logfile)

  iES_mat <- do.call(rbind, results[["scores"]])

  rownames(iES_mat) <- pathwaynames[ind]
  colnames(iES_mat) <- colnames(exp)

  return(iES_mat)
}