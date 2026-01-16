rprior<-function(n, prior, vtype="us"){

    if(!vtype%in%c("us", "idh")){
   	   stop("only 'us' and 'idh' are allowed as arguments for vtype")
    }

    if(!is.list(prior)){
    	stop("prior should be a list with elements V, nu and (optionally) alpha.mu and alpha.V")
    }

    if(!is.null(prior$fix)){
        stop("simulating with non-null fix is not yet implemented - nag me: j.hadfield@ed.ac.uk")
    }
     	
    if(is.null(prior$V)){
    	stop("V must be specified in prior")
    }

    if(is.null(prior$nu)){
    	stop("nu must be specified in prior")
    }

    if(length(prior$V)==1){
       V<-as.matrix(prior$V)
    }else{
       if(!is.matrix(prior$V)){
       	stop("V should be a matrix")
       }
    }

    if(prior$nu<=0){
      stop("improper prior: nu<=0")
    }

    if(!is.positive.definite(prior$V)){
    	stop("improper prior: V is not positive-definite")
    }
    
    is.IW<-TRUE

	if(!is.null(prior$alpha.mu)){
		is.IW<-FALSE
		if(is.null(prior$alpha.V)){
			stop("if alpha.mu is non-zero then alpha.V must be specified in prior")
	    }
	    if(!is.positive.definite(prior$alpha.V)){
	    	stop("improper prior: alpha.V is not positive-definite")
	    }
	    if(any(dim(prior$alpha.V)!=dim(prior$V))){
	    	stop("V and alpha.V should have the same dimension")
	    }
	     if(nrow(prior$alpha.V)!=length(prior$alpha.mu)){
	    	stop("alpha.mu should have the same length as the dimensions of alpha.V")
	    }
    }

    if(is.IW){
       if(vtype=="us"){	
         v<-rIW(V=prior$V, nu=prior$nu, n=n)
       }
       if(vtype=="idh"){
       	 v<-matrix(NA, n, nrow(prior$V)^2)
       	 for(k in 1:nrow(prior$V)){
       	    v[,k]<-rIW(V=as.matrix(prior$V[k,k]), nu=prior$nu, n=n)
       	 }
       }  
    }else{
       v<-matrix(NA, n, nrow(prior$V)^2)
       if(vtype=="us"){	
	       for(i in 1:n){
		       x<-MASS::mvrnorm(1, prior$alpha.mu, prior$alpha.V)
		       Vnew<-diag(x)%*%prior$V%*%diag(x)
		       v[i,]<-rIW(V=Vnew, nu=prior$nu)
	       }
	    }
	    if(vtype=="idh"){	
	       for(i in 1:n){
		       for(k in 1:nrow(prior$V)){
		       	 x<-rnorm(1, prior$alpha.mu[k], prior$alpha.V[k,k])
		         Vnew<-as.matrix(prior$V[k,k]*x^2)
		       	 v[i,k]<-rIW(V=Vnew, nu=prior$nu)
		       }
	       }
	    }

    }
    return(v)
}
