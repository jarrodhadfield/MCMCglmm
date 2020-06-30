"me"<-function(formula, error=NULL, group=NULL, type="classical"){

  
   if(!type%in%c("dberkson", "classical", "berkson", "dclassical")){stop("type must be one of 'classical', 'berkson', 'dclassical' or 'dberkson'")}
   if(type!="dberkson"){stop(paste("Sorry, meausrement error type", type, "not yet implemented"))}
   if(class(formula)!="formula"){stop("formula not passed to formula in me")}
   if(type!="dberkson" & is.null(error)){stop("Continuous covariates with mesurement error requie a standard error for 'error'")}

   formula<-update.formula(formula, ~.-1)

   X<-model.matrix(formula)
  
   if(type=="dberkson" & (any(rowSums(X)<0) | any(rowSums(X)>1) | any(X<0))){
     stop("for type='dberkson' the me formula must define probabilities of an outcome (i.e. elements and rowsums of the design matrix must lie between 0 and 1)")
   }
   if((type=="classical" | type=="berkson") & any(error<0)){
     stop("'error' must be non-negative in me terms")
   }   

   if(any(X==0)){
     X[which(X==0)]<-1e-18
   }
   if(is.null(group)){
     if(type=="dberkson"){
       attr(X, "me_prior_prob")<-cbind(1-rowSums(X), X)
     }else{
       attr(X, "me_prior_prob")<-cbind(X, error)
     }   
     attr(attr(X, "me_prior_prob"), "group")<-1:nrow(X)
   }else{
     if(length(group)!=nrow(X)){stop("group wrong length")}
     if(type=="dberkson"){
       attr(X, "me_prior_prob")<-cbind(1-rowSums(X), X)[match(unique(group), group),]
     }else{
       attr(X, "me_prior_prob")<-cbind(X, error)[match(unique(group), group),]
     } 
     attr(attr(X, "me_prior_prob"), "group")<-match(group, unique(group))
   }
   attr(X, "type")<-type
   X
}


