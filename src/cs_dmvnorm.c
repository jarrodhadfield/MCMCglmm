#include "MCMCglmm.h"

double cs_dmvnorm(const cs *beta,  const cs *mu, double ldet, const cs *Minv){

  double llik = 0.0;

  int i, j;
  double dev[Minv->m];
  double tmp;
 
  for(i=0; i<Minv->m; i++){
   dev[i] = beta->x[i]-mu->x[i];
  }

  for(i=0; i<Minv->m; i++){
    tmp= 0.0;
    for(j=0; j<Minv->m; j++){
      tmp += dev[j]*(Minv->x[j+i*Minv->m]);
    }    
    llik-= tmp*dev[i];
  }
  
  llik -= LPIx2*(Minv->m);
  llik -= ldet;
  llik /= 2.0;

  return (llik);
}                

