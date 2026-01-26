rprior<-function(n, prior, vtype="us", k=NULL){

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
       prior$V<-as.matrix(prior$V)
    }else{
       if(!is.matrix(prior$V)){
       	stop("V should be a matrix")
       }
       if(nrow(prior$V)!=ncol(prior$V)){
				stop("V should be square")
       }
    }

    if(prior$nu<=0){
      stop("improper prior: nu<=0")
    }

    if(!is.positive.definite(prior$V)){
    	stop("improper prior: V is not positive-definite")
    }
    

   if(grepl("ante", vtype)){

   	if(is.null(prior$beta.mu)){
    	  stop("beta.mu must be specified in the prior if vtype is ante")
      }
      if(is.null(prior$beta.V)){
    	  stop("beta.V must be specified in the prior if vtype is ante")
      }
      if(length(prior$beta.V)==1){
        prior$beta.V<-as.matrix(prior$beta.V)
      }else{
        if(!is.matrix(prior$beta.V)){
       	   stop("beta.V should be a matrix")
        }
        if(nrow(prior$beta.V)!=ncol(prior$beta.V)){
				stop("beta.V should be square")
        }
      }
	   if(!is.positive.definite(prior$beta.V)){
	     stop("improper prior: beta.V is not positive-definite")
	   }
      if(length(prior$beta.mu)!=nrow(prior$beta.V)){stop("beta.V should have the same dimensions as the length of beta.mu")}

   	lag<-as.numeric(gsub("[a-z]", "", vtype))
   	clag<-grepl("antec", vtype)
   	cvar<-grepl("ante.*v$", vtype)
   	vtype<-"ante"

   	if(cvar & clag){
   	   if(is.null(k)){
   		   stop("with ante models where all parameters are time homogeneous (e.g. antec1v) the size of the covariance matrix needs to be specified in the argument k")
   		}   
   	}else{
   		if(!cvar){
            k<-nrow(prior$V)
   		}else{
            k<-(length(prior$beta.mu)+sum(1:lag))/lag
   		}
   	}
   }

   if(nrow(prior$V)==1 & vtype=="us"){vtype="idh"}

   if(!vtype%in%c("us", "idh", "ante")){
   	   stop("only 'us', 'idh' and 'ante'-types are allowed as arguments for vtype")
   }
    
   is.IW<-TRUE

	if(!is.null(prior$alpha.mu)){
		is.IW<-FALSE
		if(is.null(prior$alpha.V)){
			stop("if alpha.mu is non-zero then alpha.V must be specified in prior")
	    }
	    if(length(prior$alpha.V)==1){
       	prior$alpha.V<-as.matrix(prior$alpha.V)
    	 }else{
	       if(!is.matrix(prior$alpha.V)){
	       	stop("alpha.V should be a matrix")
	       }
	       if(nrow(prior$alpha.V)!=ncol(prior$alpha.V)){
					stop("alpha.V should be square")
	       }
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
       if(vtype%in%c("idh", "ante")){
       	 v<-matrix(0, n, nrow(prior$V)^2)
       	 for(j in 1:nrow(prior$V)){
       	    v[,1+(j-1)*(nrow(prior$V)+1)]<-rIW(V=as.matrix(prior$V[j,j]), nu=prior$nu, n=n)
       	 }
       }  
    }else{
       v<-matrix(0, n, nrow(prior$V)^2)
       if(vtype=="us"){	
	       for(i in 1:n){
		       x<-MASS::mvrnorm(1, prior$alpha.mu, prior$alpha.V)
		       Vnew<-diag(x)%*%prior$V%*%diag(x)
		       v[i,]<-rIW(V=Vnew, nu=prior$nu)
	       }
	    }
	    if(vtype%in%c("idh", "ante")){	
	       for(i in 1:n){
		       for(j in 1:nrow(prior$V)){
		       	 x<-rnorm(1, prior$alpha.mu[j], prior$alpha.V[j,j])
		          Vnew<-as.matrix(prior$V[j,j]*x^2)
		       	 v[i,1+(j-1)*(nrow(prior$V)+1)]<-rIW(V=Vnew, nu=prior$nu)
		       }
	       }
	    }
    }

    if(vtype=="ante"){

      
      beta<-MASS::mvrnorm(n=n, mu=prior$beta.mu, Sigma=prior$beta.V)	
   
      if(n==1){
        beta<-matrix(beta, nrow=lag)
      }
      
      I<-diag(k)

      ante.v<-matrix(NA, n, k^2)

    	for(i in 1:n){
    		if(cvar){
         	V_epsilon<-diag(k)*v[i]
         }else{
         	V_epsilon<-matrix(v[i,], k, k)
         }
         if(clag){
         	B<-Matrix::bandSparse(k, k=-(lag:1), diagonals=sapply(lag:1, function(x){rep(beta[i,x], k-x)}))
         }else{
         	B<-Matrix::bandSparse(k, k=-(lag:1), diagonals=rev(split(beta[i,], rep(seq_along(k-1:lag), k-1:lag))))
         }
         ante.v[i,]<-c(as.matrix(solve(I-B)%*%V_epsilon%*%t(solve(I-B))))
	   }
    }
    
    return(if(vtype=="ante"){ante.v}else{v})
}
