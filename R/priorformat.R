## format an error message reporting expected and observed dimensions
dim_msg<-function(msg, expected, observed){
  sprintf("%s (expected %s, observed %s)", msg,
          paste(expected, collapse=" x "), paste(observed, collapse=" x "))
}

priorformat<-function(prior, start, nfl, meta, residual, vtype, index=NA, term=NULL){

       ## labels identifying this structure in error messages, e.g. "prior$G[[2]] (term 'animal')"
       gr<-if(residual==0){"G"}else{"R"}
       gr_index<-if(residual==0){index}else{residual}
       tlab<-if(is.null(term)){""}else{sprintf(" (term '%s')", term)}
       plab<-sprintf("prior$%s[[%s]]%s", gr, gr_index, tlab)
       slab<-sprintf("start$%s[[%s]]%s", gr, gr_index, tlab)

       if(is.null(prior)){
          prior<-list(V=diag(sum(nfl)), nu=0, fix=as.numeric(meta))
          if(grepl("ante", vtype[1])){
            nk<-as.numeric(gsub("[a-z]", "", vtype[1]))
            if(!grepl("antec", vtype[1])){
              nk<-nk*(nfl-(nk+1)/2)
            }
            prior$beta.mu<-rep(0,nk)
            prior$beta.V<-diag(nk)*1e+10
            if(grepl("ante.*v$", vtype[1])){
              prior$V<-diag(1)
            }
          }
       }

       if(inherits(prior, "prior_generator")){
         prior<-resolve_prior(prior, k=sum(nfl), vtype=vtype[1])
       }
       if(is.null(prior$V)){
         stop(paste0("V not specified for ", plab))
       }
       if(is.matrix(prior$V)==FALSE){
          prior$V<-as.matrix(prior$V)
       }
       if(!is.null(prior$alpha.mu) & residual!=0){
          stop(paste0("Parameter expanded priors not implemented for residual structures: ", plab))
       }
       if(is.null(prior$covu)){
         prior$covu<-FALSE
       }else{
         if(residual!=1){
           stop("covu is only an argument for the first residual structures")
         }else{
           if(!is.logical(prior$covu)){stop("covu should be logical")}
           if(prior$covu){
             if(grepl("ante", vtype[1]) & prior$covu){stop("covu=TRUE cannot be used with antedependence structures")}
           }
         }
       }
       if(!grepl("ante", vtype[1])){
         if(!is.null(prior$beta.mu)){
           stop("beta.mu in prior specification only possible with antedependence structures")
         }
         if(!is.null(prior$beta.V)){
           stop("beta.V in prior specification only possible with antedependence structures")
         }
       }
       if(grepl("ante.*v$", vtype[1])){
           if(!is.null(prior$alpha.V)){stop(paste0("parameter expanded priors not yet possible with antev structures - nag j.hadfield@ed.ac.uk: ", plab))}
           if(nrow(prior$V)!=1){stop("prior$V should be a scalar with antev structures")}
           prior$V<-diag(nfl)*as.numeric(prior$V)
       }
       if(is.null(prior$fix)){
         prior$fix<-0
       }
       pnames<-c("V", "n", "nu", "alpha.mu", "alpha.V", "fix", "beta.mu", "beta.V", "covu")
       if(any(names(prior)%in%pnames==FALSE)){
          paste(paste(names(prior)[which(names(prior)%in%pnames==FALSE)], sep=" "), " are not valid prior specifications for G/R-structures")
       }
       if(!prior$covu){   
         if(any(dim(prior$V)!=sum(nfl))){ # check not performed for covu - checked outside of priorformat
           if(vtype=="idvm" & nrow(prior$V)==1 & ncol(prior$V)==1){
              prior$V<-diag(prior$V[1,1], sum(nfl))
           }else{ 
              stop(dim_msg(paste0("V is the wrong dimension for ", plab), rep(sum(nfl), 2), dim(prior$V)))
           }
         }
       }
       if(is.positive.definite(prior$V)==FALSE){
         stop(paste0("V is not positive definite for ", plab))
       }
       if(is.null(prior$alpha.V)){
          prior$alpha.V<-prior$V*0
       }else{
         if(is.matrix(prior$alpha.V)==FALSE){
           prior$alpha.V<-as.matrix(prior$alpha.V)
         }
         if(any(dim(prior$alpha.V)!=dim(prior$V))){
           stop(dim_msg(paste0("alpha.V is the wrong dimension for ", plab), dim(prior$V), dim(prior$alpha.V)))
         }
         if(is.positive.definite(prior$alpha.V)==FALSE & all(prior$alpha.V==0)==FALSE){
           stop(paste0("alpha.V is not positive definite for ", plab))
         }
       } 
       if(is.null(prior$alpha.mu)){
          prior$alpha.mu<-matrix(0, nrow(prior$V), 1)
       }else{
         if(is.matrix(prior$alpha.mu)==FALSE){
           prior$alpha.mu<-matrix(prior$alpha.mu, length(prior$alpha.mu), 1)
         }
         if(length(prior$alpha.mu)!=nrow(prior$alpha.V)){
           stop(dim_msg(paste0("alpha.mu is the wrong length for ", plab), nrow(prior$alpha.V), length(prior$alpha.mu)))
         }
       }
       if(grepl("ante", vtype[1])){
         nk<-as.numeric(gsub("[a-z]", "", vtype[1]))
         if(!grepl("antec", vtype[1])){
           nk<-nk*(nfl-(nk+1)/2)
         }
         if(!is.null(prior$beta.mu)){
           if(length(prior$beta.mu)!=nk){stop(dim_msg(paste0("beta.mu is the wrong length for ", plab), nk, length(prior$beta.mu)))}
         }else{
           prior$beta.mu<-rep(0,nk)
         }
         if(!is.null(prior$beta.V)){
           if(is.matrix(prior$beta.V)==FALSE){
             prior$beta.V<-as.matrix(prior$beta.V)
           }
           if(nrow(prior$beta.V)!=ncol(prior$beta.V)){stop(sprintf("beta.V is not square for %s (observed %s)", plab, paste(dim(prior$beta.V), collapse=" x ")))}
           if(!is.positive.definite(prior$beta.V)){stop(paste0("beta.V is not positive definite for ", plab))}
           if(length(prior$beta.mu)!=nrow(prior$beta.V)){stop(dim_msg(paste0("beta.V is the wrong dimension for ", plab), rep(length(prior$beta.mu), 2), dim(prior$beta.V)))}
         }else{
           prior$beta.V<-diag(nk)*1e+10
         }
       }
       if(vtype[1]=="idvm" & nfl[1]>1){
          if(!prior$alpha.V[1,1]==0){stop(paste0("parameter expanded priors not yet possible with idvm structures - nag j.hadfield@ed.ac.uk: ", plab))}
          if(var(diag(prior$V))>0){stop("for idvm structures all diagonal elements of V should be the same")}
          if(prior$fix<2){         # either not fixed or not updated.
            prior$nu<-prior$nu/nfl # need to spread prior over all informative component.
          }else{
            prior$nu<-prior$nu/(prior$fix-1)
          }
       }
       if(prior$fix!=0){
         CM<-prior$V[prior$fix:nrow(prior$V),prior$fix:nrow(prior$V)]         
         if(prior$fix!=1){
           if(is.null(prior$nu)){
             stop(paste0("nu not specified for ", plab))
           }
           if(length(prior$nu)!=1){
             stop(paste("nu in", plab, "should be a scalar"))
           }
         }else{
           prior$nu=1
         }
       }else{
         if(is.null(prior$nu)){
           stop(paste0("nu not specified for ", plab))
         }
         if(length(prior$nu)!=1){
           stop(paste("nu in", plab, "should be a scalar"))
         }
       }
     
       if(is.null(start)){
         if(det(prior$V)<1e-8 & prior$fix==0){
           start<-prior$V+diag(nrow(prior$V))
         }else{
           start<-prior$V
         }
       }else{
         if(is.matrix(start)==FALSE){
           start<-as.matrix(start)
         }
         if(any(dim(start)!=sum(nfl)) ){
           stop(dim_msg(paste0("V is the wrong dimension for ", slab), rep(sum(nfl), 2), dim(start)))
         }
         if(is.positive.definite(start)==FALSE){
           stop(paste(slab, "is not positive definite"))
         }
       }

       pfix<-function(x,y){
         if((prior$fix>y | prior$fix<x) & prior$fix!=0){
           if(prior$fix>y){
             match(prior$fix, x:y, nomatch=0)
           }else{
             match(prior$fix, x:y, nomatch=1)
           }
         }else{
             match(prior$fix, x:y,, nomatch=0)
         }
       }

       prior<-mapply(x=cumsum(nfl)-(nfl-1), y=cumsum(nfl),  function(x,y){list(V=as.matrix(prior$V[x:y,x:y, drop=FALSE]), nu=prior$nu, fix=pfix(x,y), alpha.mu=prior$alpha.mu[x:y], alpha.V=as.matrix(prior$alpha.V[x:y,x:y, drop=FALSE]), beta.mu=prior$beta.mu, beta.V=prior$beta.V, covu=prior$covu)}, SIMPLIFY=FALSE)
       start<-mapply(x=cumsum(nfl)-(nfl-1), y=cumsum(nfl),  function(x,y){list(start=as.matrix(start[x:y,x:y, drop=FALSE]))}, SIMPLIFY=FALSE)

       return(list(prior=prior, start=start))
}

