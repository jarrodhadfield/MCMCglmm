resolve_prior <- function(prior_element, k=NULL, vtype=NULL) { 

  if(is.function(prior_element)){
    if(inherits(prior_element, "prior_generator")){
       if(is.null(k)){stop("k must be specified")}
       if(is.null(vtype)){stop("vtype must be specified")}
       prior_element(k, vtype)
    }else{
       stop("The only functions that can be passed to a prior element are the prior_generator functions IW, IG, F and stDA")
    }   
  }else{
    prior_element
  }
}