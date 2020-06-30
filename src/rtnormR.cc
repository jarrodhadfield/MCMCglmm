#include "MCMCglmmcc.h"

extern "C"{  

void rtnormR(
      int *n,
      double *mu,
      double *sd,			   
      double *lower,
      double *upper,
      double *rv
){         

  GetRNGstate();                                  

  int i;
  for(i=0; i<n[0]; i++){  
    rv[i] = rtnorm(mu[i], sd[i],lower[i], upper[i]);
  }

  PutRNGstate();

}
}

