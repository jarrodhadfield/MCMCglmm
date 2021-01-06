"pkk"<-function(prob, size){
 
    if(any(prob<0 | is.na(prob) | prob==Inf)){stop("prob must be positive")}
    if(size%%1!=0 | size<1){stop("size must be a non-zero integer")}

    prob<-prob/sum(prob)
	k<-length(prob)
    p<-0
    
    if(size>=k){

	    output<-.C("pkkR",
	      as.integer(k),
	      as.double(prob),
	      as.double(size),
	      as.double(p))

	    return(output[[4]])

	}else{

		return(0)
		
	}
}
