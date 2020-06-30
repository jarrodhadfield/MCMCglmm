#include "MCMCglmmcc.h"

extern "C"{  

void dcmvnormR(
      int *dim,
      double *betaP,
      double *muP,
      double *MP,
      int *keep,
      int *cond,			   
      int *nkeep,
      int *ncond,
      double *d
){         
 
    int i,j;
    cs *M, *mu, *beta;

    M = cs_spalloc(dim[0], dim[0], dim[0]*dim[0], true, false);
    mu = cs_spalloc(dim[0], 1, dim[0], true, false);
    beta = cs_spalloc(dim[0], 1, dim[0], true, false);

      for (i = 0 ; i < dim[0]; i++){
         M->p[i] = i*dim[0];
         mu->i[i] = i;
         mu->x[i] = muP[i];
         beta->i[i] = i;
         beta->x[i] = betaP[i];
         for (j = 0 ; j < dim[0]; j++){
           M->i[i*dim[0]+j] = j;
           M->x[i*dim[0]+j] = MP[i*dim[0]+j];
          }
      }
      mu->p[0] = 0;
      mu->p[1] = dim[0];
      beta->p[0] = 0;
      beta->p[1] = dim[0];
      M->p[dim[0]] = dim[0]*dim[0];

      d[0] = cs_dcmvnorm(beta,  mu, M, keep, nkeep[0], cond, ncond[0]);
     
    cs_spfree(M);
    cs_spfree(mu);
    cs_spfree(beta);


}
}
