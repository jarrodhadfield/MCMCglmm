"posterior.ante"<-function(x, vtype="ante1"){

    if (is.null(dim(x))[1]) {
        x <- as.matrix(x)
    }

    k<-sqrt(ncol(x))

    if (k%%1 != 0) {
        stop("cannot coerce rows of x into a matrix - the number of columns should be a triangular number")
    }

    
    lag<-as.numeric(gsub("[a-z]", "", vtype))
    
    if(!lag<sqrt(ncol(x))){stop("lag must be less than the dimensions of the matrix")}

    coerce.ante <- function(x, lag=lag) {
        k<-sqrt(length(x))
        V<-matrix(x, k, k)
        cholV<-chol(V)
        sd<-diag(diag(cholV))
        beta<-t(solve(cholV,sd))
        c(diag(sd)^2, unlist(sapply(1:lag, function(lag){-beta[(1:(k-lag)-1)*k+(1:(k-lag))+lag]})))
    }
    ante<-t(apply(x, 1, coerce.ante, lag=lag))

    colnames(ante)<-1:ncol(ante)

    if(grepl("v", vtype)){
        save_pos<-1
        colnames(ante)[1]<-paste0("iv.", sub(".*\\.", "", colnames(x)[1]))
    }else{
        save_pos<-1:k
        colnames(ante)[1:k]<-colnames(x)[seq(1, k^2, k+1)]

    }   
    if(grepl("c", vtype)){
        save_pos<-c(save_pos, k+1:lag)
        colnames(ante)[k+1:lag]<-paste0("b", 1:lag, ".", sub(".*\\.", "", colnames(x)[1]))
    }else{
        save_pos<-c(save_pos, k+1:sum(k-1:lag))
        colnames(ante)[k+1:sum(k-1:lag)]<-colnames(x)[unlist(sapply(1:lag, function(lag){seq(1+lag, k*(k-lag), k+1)}))]
    }
    ante<-ante[,save_pos]
    if (is.mcmc(x) == FALSE) {
        warning("posterior.ante expecting mcmc object")
        ante
    }
    else {
        as.mcmc(ante)
    }
}



