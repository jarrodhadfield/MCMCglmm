IW <- function(V = 1, nu = 0.002) {

  f<-function(k, vtype){

    if(!vtype%in%c("idh", "us", "idv", "idvm", "idhm")){stop(paste("prior_generator not yet implemented for", vtype, "structures"))}
    
    if(vtype=="us" | vtype=="idv" | vtype=="idvm"){
      return(list(V=diag(k)*(V*nu)/(nu+k-1), nu=nu+k-1))
    }
    if(vtype=="idh" | vtype=="idhm"){
       return(list(V=diag(k)*V, nu=nu))
    }    
  }
  class(f) <- c("prior_generator", "function")
  f
}

IG <- function(shape = 0.001, scale = 0.001) {

  f<-function(k, vtype) {

    if(!vtype%in%c("idh", "us", "idv", "idvm", "idhm")){stop(paste("prior_generator not yet implemented for", vtype, "structures"))}

  	nu<-shape*2
    V<-2*scale/nu
    if(vtype=="us" | vtype=="idv" | vtype=="idvm"){
      return(list(V=diag(k)*(V*nu)/(nu+k-1), nu=nu+k-1))
    }
    if(vtype=="idh" | vtype=="idhm"){
      return(list(V=diag(k)*V, nu=nu))
    }   
  }
  class(f) <- c("prior_generator", "function")
  f
}

F <- function(df2 = 1, scale = 1000) {

  f<-function(k, vtype) {

    if(!vtype%in%c("idh", "us", "idv", "idvm", "idhm")){stop(paste("prior_generator not yet implemented for", vtype, "structures"))}
    
    nu<-df2+k-1
    if(vtype=="us" | vtype=="idv" | vtype=="idvm"){
      return(list(V=diag(k)*df2/nu, nu=nu, alpha.mu=rep(0,k), alpha.V=diag(k)*scale))
    }
    if(vtype=="idh" | vtype=="idhm"){
      return(list(V=diag(k), nu=df2, alpha.mu=rep(0,k), alpha.V=diag(k)*scale))
    }   

  }
  class(f) <- c("prior_generator", "function")
  f
}


tSD <- function(df = 1, scale = sqrt(1000)) {
  
  f<-function(k, vtype) {
  
    if(!vtype%in%c("idh", "us", "idv", "idvm", "idhm")){stop(paste("prior_generator not yet implemented for", vtype, "structures"))}
    	
  	nu<-df+k-1
    if(vtype=="us" | vtype=="idv" | vtype=="idvm"){
      return(list(V=diag(k)*df/nu, nu=nu, alpha.mu=rep(0,k), alpha.V=diag(k)*scale^2))
    }
    if(vtype=="idh" | vtype=="idhm"){
      return(list(V=diag(k), nu=df, alpha.mu=rep(0,k), alpha.V=diag(k)*scale^2))
    }
  }
  class(f) <- c("prior_generator", "function")
  f
}

