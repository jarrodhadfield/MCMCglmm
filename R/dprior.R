dprior<-function(x, prior, sd=FALSE, log=FALSE){

    if(!is.list(prior)){
    	stop("prior should be a list with elements V, nu and (optionally) alpha.mu and alpha.V")
    }

    if(is.null(prior$V)){
    	stop("V must be specified in prior")
    }
    if(is.null(prior$nu)){
    	stop("nu must be specified in prior")
    }

    if(prior$nu<=0){
    	stop("improper prior nu<=0")
    }
    if(prior$V<=0){
    	stop("improper prior V<=0")
    }
    
    if(length(prior$V)!=1 | length(prior$nu)!=1){
    	stop("only scalar priors permitted: V should be a scalar")
    }else{
    	prior$V<-prior$V[1]
    }

    is.IW<-TRUE # Is it inverse-Wishart or parameter expanded

	if(!is.null(prior$alpha.mu)){
		is.IW<-FALSE
		if(length(prior$alpha.mu)!=1){
			stop("only scalar priors permitted: alpha.mu should be a scalar")
		}else{
           prior$alpha.mu<-prior$alpha.mu[1]
	    }		
		if(is.null(prior$alpha.V)){
			stop("if alpha.mu is non-zero then alpha.V must be specified in prior")
	    } 		
    }

	if(!is.null(prior$alpha.V)){
		if(length(prior$alpha.V)!=1){
			stop("only scalar priors permitted: alpha.V should be a scalar")
		}else{
           prior$alpha.V<-prior$alpha.V[1]
	    }			
		if(is.null(prior$alpha.mu)){
			stop("if alpha.V is non-zero then alpha.mu must be specified in prior")
	    } 		
    }

    if(is.IW){

        if(sd){
        	x<-x^2
        	jacobian<-(x^2)*0.5/sqrt(x)
        }else{
        	jacobian<-(x^2)
        }
        # note the Jacobian also includes the Jacobian for x->1/x for the inverse-gamma -> gamma

        scale<-prior$V*prior$nu

        density<-dgamma(1/x, shape=prior$nu/2, rate=scale/2, log=log)
        # Scalar inverse-Wishart (inverse-gamma) but Jacobian moved above
   
    }else{

    	if(sd){

			scale <- sqrt(prior$alpha.V*prior$V)
			jacobian<-scale

			ncp <- prior$alpha.mu/sqrt(prior$alpha.V)

            if(ncp==0){ # for half-t can just double the density

            	density<-dt(x/scale, df=prior$nu, ncp=ncp, log=log)

            	if(log){
            	 	density<-density+log(2)
            	}else{
                   	density<-density*2
            	}

            }else{

				density<-dt(x/scale, df=prior$nu, ncp=ncp)
				density<-density+dt(-x/scale, df=prior$nu, ncp=ncp)

				if(log){
					density<-log(density)
				}
			}
				
			# folded scaled non-central t
    	}else{	

			scale <- prior$alpha.V*prior$V
			jacobian<-scale

			ncp <- (prior$alpha.mu^2)/prior$alpha.V

			density<-df(x/scale, df1=1, df2=prior$nu, ncp=ncp, log=log)  
			# scaled non-central F
	    }		 
    }

    if(log){
        density<-density-log(jacobian)
    }else{
		density<-density/jacobian
    }   

    return(density)
}