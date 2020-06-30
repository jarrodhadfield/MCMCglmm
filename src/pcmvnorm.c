#include "MCMCglmm.h"

double pcmvnorm(const cs *predi, const cs *linki, const cs *G, int keep, double lower, double upper){


  double z, cpredi, cv;

  int i,
      j, 
      n = G->n,
      cnt,
      cnt2;

  if(n==1){   
    z = log(pnorm(upper, predi->x[0], sqrt(G->x[0]), TRUE,FALSE)-pnorm(lower, predi->x[0], sqrt(G->x[0]), TRUE,FALSE));
  }else{

    cs *S22, *S12, *cdev;
    csn *S22L;
    css *S22S;

    S22 = cs_spalloc (n-1, n-1, (n-1)*(n-1), 1, 0);
    S12 = cs_spalloc (n-1, 1, n-1, 1, 0);
    cdev = cs_spalloc (n-1, 1, n-1, 1, 0);

    cnt = 0;
    cnt2 = 0;
    for(i = 0; i<n; i++){ 
      if(i!=keep){
        S22->p[cnt2] = cnt2*(n-1);
        S12->x[cnt2] = G->x[i*n+keep];
        S12->i[cnt2] = cnt2;
        cdev->x[cnt2] = linki->x[i]-predi->x[i];
        cdev->i[cnt2] = cnt2;
        cnt2 ++;
      }
      for(j = 0; j<n; j++){
        if(i!=keep && j!=keep){
          S22->i[cnt] = j-1*(j>keep);
          S22->x[cnt] = G->x[i*n+j];
          cnt++;
        }
      }
    }
    S12->p[0] = 0;
    S12->p[1] = (n-1);
    cdev->p[0] = 0;
    cdev->p[1] = (n-1);
    S22->p[n-1] = (n-1)*(n-1);

    cpredi = predi->x[keep];
    cv =  G->x[keep*n+keep];

    S22S = cs_schol(1, S22);
    S22L = cs_chol(S22, S22S);

    cs_lsolve(S22L->L, S12->x);
    cs_ltsolve(S22L->L, S12->x);

    cnt=0;

    for(i = 0; i<n; i++){
      if(i!=keep){
        cpredi += cdev->x[cnt]*S12->x[cnt];
        cv  -= S12->x[cnt]*G->x[i*n+keep];
        cnt++;
      }
    }

     z = log(pnorm(upper, cpredi, sqrt(cv), TRUE,FALSE)-pnorm(lower, cpredi, sqrt(cv), TRUE,FALSE));

    cs_spfree(S22);
    cs_spfree(S12);
    cs_spfree(cdev);
    cs_nfree(S22L);
    cs_sfree(S22S);
  }
  return (z);
}                
