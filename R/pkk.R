"pkk"<-function(prob, size){
 
    prob<-prob/sum(prob)

    p<-0
    
    output<-.C("pkkR",
      as.integer(length(prob)),
      as.double(prob),
      as.double(size),
      as.double(p)
    )
    return(output[[4]])
}
