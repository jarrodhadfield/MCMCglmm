"summary.MCMCglmm"<-function(object, random=FALSE, antev=FALSE, ...){

  DIC<-object$DIC
  fixed.formula<-object$Fixed$formula
  nF<-object$Fixed$nfl
  nL<-object$Fixed$nll

  if(random){
    nF<-sum(rep(object$Random$nrl, object$Random$nfl))+nF
    if(nF!=dim(object$Sol)[2]){stop("random effects not saved and cannot be summarised")}    
  }

  solutions<-cbind(colMeans(object$Sol[,1:nF,drop=FALSE]), coda::HPDinterval(object$Sol[,1:nF,drop=FALSE]), effectiveSize(object$Sol[,1:nF,drop=FALSE]), 2*pmax(0.5/dim(object$Sol)[1], pmin(colSums(object$Sol[,1:nF,drop=FALSE]>0)/dim(object$Sol)[1], 1-colSums(object$Sol[,1:nF,drop=FALSE]>0)/dim(object$Sol)[1])))
  if(nL>0){
  solutions<-rbind(solutions, cbind(colMeans(object$Lambda), coda::HPDinterval(object$Lambda),effectiveSize(object$Lambda), 2*pmax(0.5/dim(object$Lambda)[1], pmin(colSums(object$Lambda>0)/dim(object$Lambda)[1], 1-colSums(object$Lambda>0)/dim(object$Sol)[1]))))
  }

  colnames(solutions)<-c("post.mean", "l-95% CI", "u-95% CI", "eff.samp", "pMCMC")

  random.formula=object$Random$formula
  residual.formula=object$Residual$formula

  gcomponents<-split.direct.sum(as.character(object$Random$formula)[2])
  rcomponents<-split.direct.sum(as.character(object$Residual$formula)[2])

  if(!antev){

    anteg<-grep("^ante.*\\(", gcomponents)

    if(length(anteg)>0){
       for(i in anteg){
          vtype<-gsub("\\(.*", "", gcomponents[i])
          lag<-as.numeric(gsub("[a-z]", "", vtype))

          ante_pos<-sum(object$Random$nrt[1:i])

          if(ante_pos==1){
            last<-0
          }else{
            last<-sum(object$Random$nfl[1:sum(object$Random$nfl[1:(ante_pos-1)]^2)])
          }

          k<-object$Random$nfl[ante_pos]

          post_ante<-posterior.ante(object$VCV[,last+1:(k^2)], vtype=vtype)

          object$VCV[,last+1:ncol(post_ante)]<-post_ante
          
          colnames(object$VCV)[last+1:ncol(post_ante)]<-colnames(post_ante)

          if(ncol(post_ante)!=(k^2)){
            object$VCV<-object$VCV[,-(last+(1+ncol(post_ante)):(k^2))]
            object$Random$nfl[ante_pos]<-sqrt(ncol(post_ante))
          }
       }
    }

    anter<-grep("^ante.*\\(", rcomponents)

    if(length(anter)>0){
       for(i in anter){
          vtype<-gsub("\\(.*", "", rcomponents[i])
          lag<-as.numeric(gsub("[a-z]", "", vtype))

          ante_pos<-sum(object$Residual$nrt[1:i])

          last<-sum(object$Random$nfl^2)

          if(ante_pos!=1){
            last<-last+sum(object$Residual$nfl[1:sum(object$Residual$nfl[1:(ante_pos-1)]^2)])
          }

          k<-object$Residual$nfl[ante_pos]

          post_ante<-posterior.ante(object$VCV[,last+1:(k^2)], vtype=vtype)

          object$VCV[,last+1:ncol(post_ante)]<-post_ante

          colnames(object$VCV)[last+1:ncol(post_ante)]<-colnames(post_ante)

          if(ncol(post_ante)!=(k^2)){
            object$VCV<-object$VCV[,-(last+(1+ncol(post_ante)):(k^2))]
            object$Residual$nfl[ante_pos]<-sqrt(ncol(post_ante))
          }
       }
    }
  }

  ngterms<-sum(object$Random$nfl^2)
  nrterms<-sum(object$Residual$nfl^2)

  covariances<-cbind(colMeans(object$VCV), coda::HPDinterval(object$VCV), effectiveSize(object$VCV))
  colnames(covariances)<-c("post.mean", "l-95% CI", "u-95% CI","eff.samp")
  if(ngterms>0){
   Gcovariances<-covariances[1:ngterms,,drop=FALSE]
  }else{
   Gcovariances<-NULL
  }
  Rcovariances<-covariances[ngterms+1:nrterms,,drop=FALSE]
  cstats<-attr(object$VCV, "mcpar")
  cstats[4]<-dim(object$VCV)[1]
  if(is.null(object$CP)){
     cutpoints<-NULL
  }else{
    cutpoints<-cbind(colMeans(object$CP), coda::HPDinterval(object$CP), effectiveSize(object$CP))
    colnames(cutpoints)<-c("post.mean", "l-95% CI", "u-95% CI", "eff.samp")

  }
  if(is.null(object$ThetaS)){
     theta_scale<-NULL
  }else{
    theta_scale<-cbind(colMeans(object$ThetaS), coda::HPDinterval(object$ThetaS), effectiveSize(object$ThetaS), 2*pmax(0.5/dim(object$ThetaS)[1], pmin(colSums(object$ThetaS>0)/dim(object$ThetaS)[1], 1-colSums(object$ThetaS>0)/dim(object$ThetaS)[1])))

      colnames(theta_scale)<-c("post.mean", "l-95% CI", "u-95% CI", "eff.samp", "pMCMC")
  }
  if(is.null(object$Random$nrt)){
    Gterms<-NULL
  }else{
    Gterms<-rep(rep(1:length(object$Random$nrt), object$Random$nrt), object$Random$nfl^2) 
  }
  Rterms<-rep(rep(1:length(object$Residual$nrt), object$Residual$nrt), object$Residual$nfl^2) 

  output<-list(DIC=DIC, fixed.formula=fixed.formula, random.formula=random.formula,residual.formula=residual.formula, solutions=solutions, Gcovariances=Gcovariances, Gterms=Gterms,  Rcovariances=Rcovariances,Rterms=Rterms, cstats=cstats,cutpoints=cutpoints, theta_scale=theta_scale)
  attr(output, "class")<-c("summary.MCMCglmm", "list")
  output
}

"print.summary.MCMCglmm"<-function (x, digits = max(3, getOption("digits") - 3), has.Pvalue=TRUE, eps.Pvalue = 1/(x$cstats[4]-1), cstats=TRUE,  ...) 
{

 if(cstats){
   cat("\n Iterations =", paste(x$cstats[1], ":", x$cstats[2], sep=""))
   cat("\n Thinning interval  =" , x$cstats[3]) 
   cat("\n Sample size  =" , x$cstats[4], "\n") 
 }

 cat("\n DIC:", x$DIC, "\n")
 if(is.null(x$random.formula)==FALSE){
   rcomponents<-split.direct.sum(as.character(x$random.formula)[2])
   for(i in 1:length(rcomponents)){
     if(i==1){
     cat(paste("\n G-structure:  ~", rcomponents[i], "\n\n", sep=""))
     }else{
     cat(paste("\n               ~", rcomponents[i], "\n\n", sep=""))
     }
     if(i%in%x$Gterms){
       print(as.data.frame(x$Gcovariance[x$Gterms==i,,drop=FALSE]), digits=digits, ...)
     }else{
       cat(" G-R structure below\n")
     }
   }
 }
 rcomponents<-split.direct.sum(as.character(x$residual.formula)[2])
 for(i in 1:length(rcomponents)){
   if(i==1){
     cat(paste("\n R-structure:  ~", rcomponents[i], "\n\n", sep=""))
   }else{
     cat(paste("\n               ~", rcomponents[i], "\n\n", sep=""))
   }
   print(as.data.frame(x$Rcovariance[x$Rterms==i,,drop=FALSE]), digits=digits, ...)
 }

 cat("\n Location effects:", paste(as.expression(x$fixed.formula)), "\n\n")
 printCoefmat(as.data.frame(x$solutions), has.Pvalue=has.Pvalue, digits=digits, eps.Pvalue=eps.Pvalue, ...)

 if(!is.null(x$cutpoints)){
   cat("\n Cutpoints:", "\n\n")
   print(as.data.frame(x$cutpoints), digits=digits, ...)
 }

  if(!is.null(x$theta_scale)){
   cat("\n Theta scale parameter:", "\n\n")
   printCoefmat(as.data.frame(x$theta_scale), has.Pvalue=has.Pvalue, digits=digits, eps.Pvalue=eps.Pvalue, ...)
 }

}

