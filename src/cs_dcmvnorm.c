#include "MCMCglmm.h"

double cs_dcmvnorm(const cs *beta,  const cs *mu, const cs *M, int *keep, int nkeep, int *cond, int ncond){

  double llik = 0.0;

  int i,
      j, 
      cnt;
  double ldet;

  cs *S11, *S22, *invS22, *invS11, *S12, *S21, *muC, *dev, *cdev;

  if(ncond>0){

    S22 = cs_spalloc (ncond, ncond, ncond*ncond, 1, 0);
    invS22 = cs_spalloc (ncond, ncond, ncond*ncond, 1, 0);
    S12 = cs_spalloc (nkeep, ncond, nkeep*ncond, 1, 0);
    cdev = cs_spalloc (ncond, 1, ncond, 1, 0);

    cnt = 0;
    for(i = 0; i<ncond; i++){ 
      S22->p[i] = i*ncond;
      for(j = 0; j<ncond; j++){
        S22->i[cnt] = j;
        S22->x[cnt] = M->x[M->p[cond[i]]+cond[j]];
        cnt++;
      }
    }
    S22->p[ncond] = ncond*ncond;

    cnt = 0;
    for(i = 0; i<ncond; i++){
      invS22->p[i] = i*ncond;
      cdev->i[i] = i;
      cdev->x[i] = beta->x[cond[i]]-mu->x[cond[i]];
      for(j = 0; j<ncond; j++){
        invS22->i[cnt] = j;
        invS22->x[cnt] = 1.0;
        cnt++;
      }
    }
    cdev->p[0] = 0;
    cdev->p[1] = ncond;
    invS22->p[ncond] = ncond*ncond;

    cs_invR(S22, invS22);

    cnt = 0;
    for(i = 0; i<ncond; i++){
      S12->p[i] = i*nkeep;
      for(j = 0; j<nkeep; j++){
        S12->i[cnt] = j;
        S12->x[cnt] = M->x[M->p[cond[i]]+keep[j]]; 
        cnt ++;
      }
    }
    S12->p[ncond] = ncond*nkeep;

    muC = cs_multiply(S12, invS22);
    S21 = cs_transpose(S12, TRUE);
    S11 = cs_multiply(muC, S21);
    dev = cs_multiply(muC, cdev);

    cnt = 0;

    for(i = 0; i<nkeep; i++){
      for(j = 0; j<nkeep; j++){
        S11->x[cnt] = M->x[M->p[keep[i]]+keep[j]]-S11->x[cnt];
        cnt ++;
      }
    }

  }else{

    S11 = cs_spalloc (nkeep, nkeep, nkeep*nkeep, 1, 0);
    dev = cs_spalloc (nkeep, 1, nkeep, 1, 0);

    cnt = 0;
    for(i = 0; i<nkeep; i++){ 
      S11->p[i] = i*nkeep;
      dev->i[i] = i;
      dev->x[i] = 0.0;
      for(j = 0; j<nkeep; j++){
        S11->i[cnt] = j;
        S11->x[cnt] = M->x[M->p[keep[i]]+keep[j]];
        cnt++;
      }
    }
    dev->p[0] = 0;
    dev->p[1] = nkeep;
    S11->p[nkeep] = nkeep*nkeep;
  }

  invS11 = cs_spalloc (nkeep, nkeep, nkeep*nkeep, 1, 0);

  cnt = 0;
  for(i = 0; i<nkeep; i++){
    invS11->p[i] = i*nkeep;
    for(j = 0; j<nkeep; j++){
      invS11->i[cnt] = j;
      invS11->x[cnt] = 1.0;
      cnt++;
    }
  }
  invS11->p[nkeep] = nkeep*nkeep;

  ldet = log(cs_invR(S11, invS11));
 
  for(i = 0; i<nkeep; i++){
    dev->x[i]  = beta->x[keep[i]]-dev->x[i]-mu->x[keep[i]];
  }

  for(i=0; i<nkeep; i++){
    for(j=0; j<nkeep; j++){
      llik -= dev->x[j]*(invS11->x[j+i*nkeep])*dev->x[i];
    }    
  }
 
  llik -= LPIx2*(nkeep);
  llik -= ldet;
  llik /= 2.0;

  cs_spfree(S11);
  cs_spfree(invS11);
  cs_spfree(dev);

  if(ncond>0){
    cs_spfree(S22);
    cs_spfree(invS22);
    cs_spfree(S12);
    cs_spfree(S21);
    cs_spfree(muC);
    cs_spfree(cdev);
  }

  return (llik);
}                

