"predict.MCMCglmm"<-function(object, newdata=NULL, marginal=object$Random$formula, type="response", interval="none", level=0.95, it=NULL, posterior="all", verbose=FALSE, approx="numerical", ...){

  rm.obs<-c()

  if(interval=="prediction"){
      post.pred<-t(simulate(object=object, nsim=nrow(object$Sol), newdata=newdata, marginal=marginal, type=type, it=it, posterior=posterior, verbose=verbose))

  }else{

    if(type%in%c("response", "terms")==FALSE){stop("type must be response or terms")}
    if(interval%in%c("none", "confidence", "prediction")==FALSE){stop("interval must be none, confidence or prediction")}

    if(!is.null(posterior)){
      if(!posterior%in%c("distribution", "mean", "mode", "all")){
        stop("posterior argument must be either distribution, mean, mode or all")
      }
    }

    if(!is.null(marginal)){
      if(!is(marginal, "formula")){stop("marginal should be NULL or a formula")}
    }

    if(!is.null(it)){
      if(length(it)>1){stop("it should be an integer")}
      if(it>nrow(object$Sol) | it<1){stop("it should be less than or equal to the number of iterations")}
    }

    rcomponents<-split.direct.sum(as.character(object$Random$formula)[2])
    mcomponents<-split.direct.sum(as.character(marginal)[2])
    
    if(is.null(object$meta)){object$meta<-FALSE}
      
    if((length(rcomponents)+object$meta)!=length(object$Random$nrt)){stop("if mev was used, add my_model$meta=TRUE and rerun. If not, the model is a covu model and predictions are not yet implemented for this type of model")}

    if(any(mcomponents%in%rcomponents==FALSE)){stop("marginal formula does not correspond to model formula")}

    marginalise<-rep(as.numeric(rcomponents%in%mcomponents), object$Random$nrt)

    if(!is.null(newdata)){
      suppressWarnings(object2<-MCMCglmm(fixed=object$Fixed$formula, random=object$Random$formula, rcov=object$Residual$formula, family=object$Residual$original.family, data=newdata, nitt=1, thin=1, burnin=0, ginverse=object$ginverse, verbose=FALSE, pr=any(marginalise==0)))
      find.fixed<-match(colnames(object2$Sol)[1:object2$Fixed$nfl], colnames(object$Sol))
      find.random<-match(colnames(object2$Sol)[-c(1:object2$Fixed$nfl)], colnames(object$Sol))
      if(any(is.na(find.fixed))){stop("model for newdata has fixed effects not present in original model")}
      if(any(!colnames(object$Sol)[1:object$Fixed$nfl]%in%colnames(object2$Sol)[1:object2$Fixed$nfl])){stop("model for newdata has fixed effects absent from the original model")}
      # Jarrod; previously this error catching was not implemented, but a waraning was issued if verbose=TRUE (L53).
      if(verbose){
        if(any(is.na(find.random))){
          missing.random<-colnames(object2$Sol)[which(is.na(find.random))+object2$Fixed$nfl]
          warning(paste("model for newdata has random effects not present in original model:", paste(missing.random, collapse=", ")))
        }
        missing.fixed<-which(!colnames(object$Sol)[1:object$Fixed$nfl]%in%colnames(object2$Sol)[1:object2$Fixed$nfl])
        if(length(missing.fixed)>0){
          missing.fixed<-colnames(object$Sol)[1:object$Fixed$nfl][missing.fixed]
          warning(paste("original model has fixed effects not present in newdata:", paste(missing.fixed, collapse=", ")))
        }
      }
      object2$Sol<-object$Sol[,c(find.fixed, find.random),drop=FALSE]
      find.vcv<-match(colnames(object2$VCV), colnames(object$VCV))
      if(any(is.na(find.vcv))){stop("model for newdata has (co)variance terms not in original model")}
      object2$VCV<-object$VCV[,find.vcv,drop=FALSE]
      if(!is.null(object2$CP)){
        find.cp<-match(colnames(object2$CP), colnames(object$CP))
        if(any(is.na(find.cp))){stop("model for newdata has cutpoints not in original model")}
        object2$CP<-object$CP[,find.cp,drop=FALSE]
      }
      object<-object2
      rm(object2)
    }

    if(posterior=="mean"){
      object$VCV<-matrix(colMeans(object$VCV), 1, ncol(object$VCV))
	    object$Sol<-matrix(colMeans(object$Sol), 1, ncol(object$Sol))
      if(!is.null(object$CP)){
        object$CP<-matrix(colMeans(object$CP), 1, ncol(object$CP))
      }
      if(!is.null(object$Lambda)){
        object$Lambda<-matrix(colMeans(object$Lambda), 1, ncol(object$Lambda))
      }
      it<-1
    }
    if(posterior=="mode"){
      object$VCV<-matrix(posterior.mode(object$VCV, ...), 1, ncol(object$VCV))
	    object$Sol<-matrix(posterior.mode(object$Sol, ...), 1, ncol(object$Sol))
      if(!is.null(object$CP)){
        object$CP<-matrix(posterior.mode(object$CP, ...), 1, ncol(object$CP))
      }
      if(!is.null(object$Lambda)){
        object$Lambda<-matrix(posterior.mode(object$Lambda, ...), 1, ncol(object$Lambda))
      }
      it<-1
    }
    if(is.null(it)){
      if(posterior=="distribution"){
        it<-sample(nrow(object$Sol), 1)
      }
      if(posterior=="all"){
        it<-1:nrow(object$Sol)
      }
    }
    object$Sol<-object$Sol[it,,drop=FALSE]
    object$VCV<-object$VCV[it,,drop=FALSE]
    if(!is.null(object$Lambda)){
      object$Lambda<-object$Lambda[it,,drop=FALSE]
    }
    if(!is.null(object$CP)){
      object$CP<-object$CP[it,,drop=FALSE]
    }

    if(is.null(object$X)){
      stop("fixed effect design matrix not saved: pass saveX=TRUE to MCMCglmm")
    }

    if(any(marginalise==0) & dim(object$Sol)[2]==dim(object$X)[2]){
      stop("posterior distribution of random effects not saved: pass pr=TRUE to MCMCglmm")
    }
    if(any(marginalise==0) & is.null(object$Z)){
      stop("random effect design matrix not saved: pass saveZ=TRUE to MCMCglmm")
    }


    if(is.null(object$Random$nfl)==FALSE){  # there are random effects
      st<-c(1,cumsum(object$Random$nrl*object$Random$nfl)+1)  # starting column for random effects of each component 
      st<-st[-length(st)]
      end<-cumsum(object$Random$nrl*object$Random$nfl)        # ending column for random effects of each component 
      keep<-unlist(mapply(st[which(marginalise==0)], end[which(marginalise==0)], FUN=":"))    # random effects to be kept
    }else{
      keep<-NULL
    }

    if(!is.null(newdata) & !is.null(object$Random$nfl)){
      missing.random<-which(is.na(object$Sol[1,]))
      if(posterior=="mean" | posterior=="mode"){
        if(any(is.na(object$Sol))){
          object$Sol[missing.random]<-0
        }
      }else{
        dv<-1:ncol(object$Sol)
        cnt<-0
        for(j in 1:length(object$Random$nfl)){
          nfl<-object$Random$nfl[j]
          nrl<-object$Random$nrl[j]
          dv[object$Fixed$nfl+cnt+1:(nfl*nrl)]<-rep(diag(matrix(1:(nfl^2),nfl,nfl)), each=nrl)
          cnt<-cnt+(nfl*nrl)
        }
        if(any(is.na(object$Sol))){
          object$Sol[,missing.random,drop=FALSE]<-rnorm(nrow(object$Sol)*length(missing.random), 0, sqrt(object$VCV[,dv[missing.random]]))
        }
      }
    }

    object$Sol<-object$Sol[,c(1:object$Fixed$nfl, object$Fixed$nfl+keep),drop=FALSE]

    W<-cbind(object$X, object$Z)
    W<-W[,c(1:object$Fixed$nfl, object$Fixed$nfl+keep), drop=FALSE]
  
    if(nrow(object$Sol)==1){
       post.pred<-t(as.matrix(W%*%t(object$Sol)))
    }else{
       post.pred<-t(apply(object$Sol,1, function(x){(W%*%x)@x}))
    }
    if(type=="response"){
 
      if(any(object$family!="gaussian" & object$family!="cengaussian" & object$family!="ncst")){
         post.var<-buildV(object, marginal=marginal, diag=TRUE, it=NULL, posterior="all", ...)
      }

      ordering<-object$ZR@i+1

      post.pred<-post.pred[,ordering, drop=FALSE]         
      post.var<-post.var[,ordering, drop=FALSE]
      object$family<-object$family[ordering]
      object$y.additional<-object$y.additional[ordering,]
      object$error.term<-object$error.term[ordering]
      # when predicting on the response scale best to order everything in terms of R-order 

      normal.evd<-function(mu, v, approx){
        if(approx=="numerical" || approx=="taylor2"){
          if(approx=="numerical"){ 
            int.foo<-function(x, mu, v){exp(-exp(x))*dnorm(x, mu, sqrt(v))}
            return(integrate(int.foo, qnorm(1e-6, mu,sqrt(v)), qnorm(1-1e-6, mu,sqrt(v)), mu,v)[[1]])
          }else{
            return(exp(-exp(mu))+0.5*v*exp(-exp(mu)+mu)*(exp(mu)-1))
          }
        }else{  
          stop(paste(approx, "approximation not implemented for this response"))
        }  
      }
      normal.ztp<-function(mu, v, approx){
        if(approx=="numerical" || approx=="taylor2"){
          if(approx=="numerical"){
            int.foo<-function(x, mu, v){exp(x)/(1-exp(-exp(x)))*dnorm(x, mu, sqrt(v))}
            return(integrate(int.foo, qnorm(0.0001, mu,sqrt(v)), qnorm(0.9999, mu,sqrt(v)), mu,v)[[1]])
          }else{
            return(exp(mu)/(1-exp(-exp(mu)))+0.5*v*exp(exp(mu)+mu)*(1-2*exp(exp(mu))+exp(2*exp(mu))+3*exp(mu)+2*exp(2*mu)-3*exp(exp(mu)+mu)+exp(exp(mu)+2*mu))/(exp(exp(mu))-1)^3)
          } 
        }else{  
          stop(paste(approx, "approximation not implemented for this response"))
        }  
      }  
      normal.ztb<-function(mu, v, size, approx){
        if(approx=="numerical" || approx=="taylor2"){
          if(approx=="numerical"){
            int.foo<-function(x, mu, v){(size*plogis(x)/(1-(1-plogis(x))^size))*dnorm(x, mu, sqrt(v))}
            return(integrate(int.foo, qnorm(0.0001, mu,sqrt(v)), qnorm(0.9999, mu,sqrt(v)), mu,v)[[1]])
          }else{
            stop(paste(approx, "taylor approximation not implemented for this response"))
          } 
        }else{  
          stop(paste(approx, "approximation not implemented for this response"))
        }  
      }  
      normal.logistic<-function(mu, v, approx){
        if(approx=="numerical"){
          int.foo<-function(x, mu, v){plogis(x)*dnorm(x, mu, sqrt(v))}
          return(integrate(int.foo, qnorm(0.0001, mu,sqrt(v)), qnorm(0.9999, mu,sqrt(v)), mu,v)[[1]])
        }
        if(approx=="diggle"){
          return(plogis(mu/sqrt(1+v*(16*sqrt(3)/(15*pi))^2)))
        }
        if(approx=="mcculloch"){
          return(plogis(mu-0.5*v*tanh(mu*(1+2*exp(-0.5*v))/6)))
        }
        if(approx=="taylor2"){
          return(plogis(mu)-0.5*v*exp(mu)*(exp(mu)-1)/(1+exp(mu))^3)
        }
      }
  
      normal.clogistic<-function(mu, v, size, approx){
        if(approx=="numerical" || approx=="taylor2"){
          if(approx=="numerical"){
            int.foo<-function(x, mu, v, size){(1-((1-plogis(mu))^size))*dnorm(x, mu, sqrt(v))}
            return(integrate(int.foo, qnorm(0.0001, mu,sqrt(v)), qnorm(0.9999, mu,sqrt(v)), mu,v, size)[[1]])
          }else{
            return((1-((1-plogis(mu))^size))-0.5*v*exp(mu)*size*(size*exp(mu)-1)/(1+exp(mu))^(2+size))
          }
        }else{
          stop(paste(approx, "approximation not implemented for this response"))
        }  
      }
      
      normal.multilogistic<-function(mu, v, approx){
        if(approx=="numerical"){
          int.foo<-function(x, mu, v, i){(exp(x[i])/(1+sum(exp(x))))*prod(dnorm(x, mu, sqrt(v)))}
          res<-1:length(mu)
          for(i in 1:length(mu)){
            res[i]<-cubature::adaptIntegrate(int.foo, qnorm(0.0001, mu,sqrt(v)), qnorm(0.9999, mu,sqrt(v)), mu=mu, v=v, i=i)[[1]]
          }
          return(res)
        }else{
          stop(paste(approx, "approximation not implemented for this response"))
        }
      } 

      if(!is.null(object$Lambda)){
        if(any(object$family!="gaussian")){stop("sorry - prediction not yet available for sir models applied to non-Gaussian data")}
        # bear in mind that if this is implemented the post.var is quadratic in the structural parameters
        I<-Diagonal(nrow(object$XL))
        post.pred<-t(sapply(1:nrow(object$Sol), function(i){as.vector(solve(I-object$XL%*%kronecker(object$Lambda[i,], I), post.pred[i,]))}))
      }
      if(!is.null(object$ThetaS)){stop("sorry - prediction not yet available for theta_scale models")}

      if(any(object$family%in%c("poisson","cenpoisson"))){
        keep<-which(object$family%in%c("poisson","cenpoisson"))
        post.pred[,keep]<-exp(post.pred[,keep]+0.5*post.var[,keep])
      }

      if(any(object$family%in%c("ztpoisson"))){
        keep<-which(object$family%in%c("ztpoisson"))
        post.pred[,keep]<-mapply(post.pred[,keep], post.var[,keep], FUN=function(mu,v){normal.ztp(mu,v, approx)})
      }

      if(any(object$family%in%c("exponential","cenexponential","geometric","cengeometric"))){
        keep<-which(object$family%in%c("exponential","cenexponential","geometric","cengeometric"))
        post.pred[,keep]<-exp(-post.pred[,keep]+0.5*post.var[,keep])
      }

      if(any(object$family%in%c("ordinal", "threshold"))){      

        nord<-unique(object$error.term[which(object$family%in%c("ordinal", "threshold"))])
        cp.names<-substr(colnames(object$CP), 10, regexpr("\\.[1-9]$|\\.[1-9][0-9]$", colnames(object$CP))-1)
        cp.names<-match(cp.names,unique(cp.names))

        for(k in 1:length(nord)){

          keep<-which(object$family%in%c("ordinal", "threshold") & object$error.term==nord[k])

          which.ord<-object$Residual$mfac[nord[k]]+1 

          CP<-cbind(-Inf, 0, object$CP[,which(cp.names==which.ord), drop=FALSE], Inf)
          q<-matrix(0, dim(post.pred)[1], length(keep))

          is.ordinal<-as.numeric(object$family[keep[1]]=="ordinal")

          for(i in 2:(dim(CP)[2]-1)){
            q<-q+(pnorm(CP[,i+1]-post.pred[,keep],0,sqrt(post.var[,keep]+is.ordinal))-pnorm(CP[,i]-post.pred[,keep],0,sqrt(post.var[,keep]+is.ordinal)))*(i-1)
          }

          post.pred[,keep]<-q
          rm(q)
        } 
      }

      if(any(grepl("nzbinom", object$family))){
        keep<-grep("nzbinom", object$family)
        size<-t(matrix(object$y.additional[keep,1], length(keep),nrow(post.pred)))
        post.pred[,keep]<-mapply(post.pred[,keep], post.var[,keep], size, FUN=function(mu,v,size){normal.clogistic(mu,v,size, approx)})
      }
      
      if(any(grepl("ncst", object$family))){
        keep<-grep("ncst", object$family)
        scale<-t(matrix(object$y.additional[keep,1], length(keep),nrow(post.pred)))
        df<-t(matrix(object$y.additional[keep,2], length(keep),nrow(post.pred)))
        post.pred[,keep]<-exp(0.5*log(df/2)+lgamma(df/2-1/2)-lgamma(df/2))*post.pred[,keep]
      }


      if(any(grepl("multinomial", object$family))){

        keep<-grep("multinomial", object$family)
        eterms<-unique(object$error.term[keep])
        ncat<-rle(object$Residual$mfac[eterms])
        ncat<-rep(ncat$values + 2, ncat$lengths/(ncat$values+1))
        # recover the number of categories of each multinomial supertrait from mfac

        et<-0
        for(i in 1:length(ncat)){
          keep<-which(grepl("multinomial", object$family) & object$error.term%in%eterms[et+1:(ncat[i]-1)])
          size<-object$y.additional[keep[seq(1, by=ncat[i]-1, length(keep))],1]
       
          for(j in 1:nrow(post.pred)){
            prob<-matrix(post.pred[j,keep], ncol=ncat[i]-1)
            pvar<-matrix(post.var[j,keep], ncol=ncat[i]-1)
            if(ncat[i]==1){
              post.pred[j,keep]<-t(sapply(1:nrow(prob), function(x){normal.logistic(prob[x,], pvar[x,], approx)}))*size
            }else{
              post.pred[j,keep]<-t(sapply(1:nrow(prob), function(x){normal.multilogistic(prob[x,], pvar[x,], approx)}))*size
            }
          }
          et<-cumsum(ncat[1:i]-1) # end point of the multinomial super.trait in eterms
        }
      }


      if(any(grepl("hupoisson", object$family))){
        keep<-grep("hupoisson", object$family)
        eterms<-unique(object$error.term[keep])
        eterms<-eterms[seq(1, by=2, length(eterms))]
        for(i in 1:length(eterms)){
          keep1<-which(grepl("hupoisson", object$family) & object$error.term%in%(eterms[i]+1))
          post.pred[,keep1]<-mapply(post.pred[,keep1], post.var[,keep1], FUN=function(mu,v){1-normal.logistic(mu,v, approx)})
          keep2<-which(grepl("hupoisson", object$family) & object$error.term%in%(eterms[i]))
          post.pred[,keep2]<-post.pred[,keep1]*mapply(post.pred[,keep2], post.var[,keep2], FUN=function(mu,v){normal.ztp(mu,v, approx)})
          rm.obs<-c(rm.obs, keep1)
        }  
      }

      if(any(grepl("zapoisson", object$family))){
        keep<-grep("zapoisson", object$family)
        eterms<-unique(object$error.term[keep])
        eterms<-eterms[seq(1, by=2, length(eterms))]
        for(i in 1:length(eterms)){
          keep1<-which(grepl("zapoisson", object$family) & object$error.term%in%(eterms[i]+1))
          post.pred[,keep1]<-mapply(post.pred[,keep1], post.var[,keep1], FUN=function(mu,v){1-normal.evd(mu,v, approx)})
          keep2<-which(grepl("zapoisson", object$family) & object$error.term%in%(eterms[i]))
          post.pred[,keep2]<-post.pred[,keep1]*mapply(post.pred[,keep2], post.var[,keep2], FUN=function(mu,v){normal.ztp(mu,v, approx)})
          rm.obs<-c(rm.obs, keep1)
        }  
      }

      if(any(grepl("zipoisson", object$family))){
        keep<-grep("zipoisson", object$family)
        eterms<-unique(object$error.term[keep])
        eterms<-eterms[seq(1, by=2, length(eterms))]
        for(i in 1:length(eterms)){
          keep1<-which(grepl("zipoisson", object$family) & object$error.term%in%(eterms[i]+1))
          post.pred[,keep1]<-mapply(post.pred[,keep1], post.var[,keep1], FUN=function(mu,v){1-normal.logistic(mu,v, approx)})
          keep2<-which(grepl("zipoisson", object$family) & object$error.term%in%(eterms[i]))
          post.pred[,keep2]<-post.pred[,keep1]*exp(post.pred[,keep2]+post.var[,keep2]/2)
          rm.obs<-c(rm.obs, keep1)
        }  
      }

      if(any(grepl("hubinomial", object$family))){
        keep<-grep("hubinomial", object$family)
        eterms<-unique(object$error.term[keep])
        eterms<-eterms[seq(1, by=2, length(eterms))]
        for(i in 1:length(eterms)){
          keep1<-which(grepl("hubinomial", object$family) & object$error.term%in%(eterms[i]+1))
          post.pred[,keep1]<-mapply(post.pred[,keep1], post.var[,keep1], FUN=function(mu,v){1-normal.logistic(mu,v, approx)})
          keep2<-which(grepl("hubinomial", object$family) & object$error.term%in%(eterms[i]))
          size<-object$y.additional[keep2,1]
          post.pred[,keep2]<-post.pred[,keep1]*mapply(post.pred[,keep2], post.var[,keep2], rep(size, each=nrow(post.pred)), FUN=function(mu,v, size){normal.ztb(mu,v, size, approx)})
          rm.obs<-c(rm.obs, keep1)
        }  
      }

      if(any(grepl("zibinomial", object$family))){
        keep<-grep("zibinomial", object$family)
        eterms<-unique(object$error.term[keep])
        eterms<-eterms[seq(1, by=2, length(eterms))]
        for(i in 1:length(eterms)){
          keep1<-which(grepl("zibinomial", object$family) & object$error.term%in%(eterms[i]+1))
          post.pred[,keep1]<-mapply(post.pred[,keep1], post.var[,keep1], FUN=function(mu,v){1-normal.logistic(mu,v, approx)})
          keep2<-which(grepl("zibinomial", object$family) & object$error.term%in%(eterms[i]))
          size<-object$y.additional[keep2,1]
          post.pred[,keep2]<-size*post.pred[,keep1]*mapply(post.pred[,keep2], post.var[,keep2], FUN=function(mu,v){normal.logistic(mu,v, approx)})
          rm.obs<-c(rm.obs, keep1)
        }  
      }

      post.pred<-post.pred[,order(ordering), drop=FALSE]
      rm.obs<-match(rm.obs, (1:length(ordering))[order(ordering)])
      # put back in original order

    }
  }
  
  if(length(rm.obs)>0){
    post.pred<-post.pred[,-rm.obs, drop=FALSE]
  }

  if(nrow(post.pred)>1){
    pred<-matrix(colMeans(post.pred), dim(post.pred)[2],1)
  }else{
    pred<-t(post.pred)
  }  

  if(interval!="none"){
    pred<-cbind(pred, coda::HPDinterval(mcmc(post.pred), prob=level))   
    colnames(pred)<-c("fit", "lwr", "upr")
  }
  
  rownames(pred)<-1:dim(pred)[1]

  return(pred)
}  
