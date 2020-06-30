#include "MCMCglmmcc.h"

extern "C"{  

void rtcmvnormR(
      int *n,
      double *muP,
      double *mu2P,
      double *GP,
      int *keep,
      int *dim,			   
      double *lower,
      double *upper,
      double *rv
){         


    int i,j;
    cs *G, *mu, *mu2;

    G = cs_spalloc(dim[0], dim[0], dim[0]*dim[0], true, false);
    mu = cs_spalloc(dim[0], 1, dim[0], true, false);
    mu2 = cs_spalloc(dim[0], 1, dim[0], true, false);

    for (i = 0 ; i < dim[0]; i++){
       G->p[i] = i*dim[0];
       mu->i[i] = i;
       mu->x[i] = muP[i];
       mu2->i[i] = i;
       mu2->x[i] = mu2P[i];
       for (j = 0 ; j < dim[0]; j++){
         G->i[i*dim[0]+j] = j;
         G->x[i*dim[0]+j] = GP[i*dim[0]+j];
       }
    }
    mu->p[0] = 0;
    mu->p[1] = dim[0];
    mu2->p[0] = 0;
    mu2->p[1] = dim[0];
    G->p[dim[0]] = dim[0]*dim[0];

    GetRNGstate();                                  

    for(i=0; i<n[0]; i++){  
      rv[i] = rtcmvnorm(mu, mu2, G, keep[0], lower[0], upper[0]);
    }

    PutRNGstate();

    cs_spfree(G);
    cs_spfree(mu);
    cs_spfree(mu2);
}
}
