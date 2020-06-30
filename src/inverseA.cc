#include "MCMCglmmcc.h"

extern "C"{  

/**************************/
/* Meuwissen and Luo 1992 */
/**************************/

void inverseA(
        int *id,      
        int *dam,       
        int *sire,         
        double *f,     
        double *dii,    
        int *iAP,              
	int *pAP,	         
	double *xAP,
        int *nAP,
	int *nzmaxAP,
	int *iTinvP,              
	int *pTinvP,	         
	double *xTinvP,
        int *nTinvP,
	int *nzmaxTinvP
){         

  int     i, j, k, cnt, sj, dj;
  double  ai;
  double  *AN = new double[2*nAP[0]];
 // needs to hold all ancestors of an individual (btu ancestors can be counted multiple times if inbreeding occurs) 
 // 2* the total number of individuals is probably safe (that observed under a pedigree formed of one full-sib mating lineage).
  double  *li = new double[nAP[0]];
  cs *Tinv, *D, *tTinv, *tTD, *A, *tA;

  for(i=0; i<nAP[0]; i++){
     li[i]=0.0;               // set l to zero
  }
  for(i=0; i<nAP[0]; i++){
     AN[i]=-1;               // set AN to zero
  }

  Tinv = cs_spalloc(nTinvP[0], nTinvP[0], nzmaxTinvP[0], true, false);  

         for (i = 0 ; i < nzmaxTinvP[0] ; i++){
           Tinv->i[i] = iTinvP[i];
           Tinv->x[i] = xTinvP[i];
         }
         for (i = 0 ; i <= nTinvP[0]; i++){
           Tinv->p[i] = pTinvP[i];
         }

  tTinv = cs_transpose(Tinv, true);

  D = cs_spalloc(nTinvP[0], nTinvP[0], nzmaxTinvP[0], true, false);  

        for (i = 0 ; i < nTinvP[0] ; i++){
           D->i[i] = i;
           D->x[i] = 1.0;
           D->p[i] = i;
         }
         D->p[nTinvP[0]] = nTinvP[0];
    
  for(i=0; i<nAP[0]; i++){  // iterate through each row of l 

    li[i] = 1.0;                   // set l_ii to one
    ai=0.0;                        // set a_ii to zero

    if(dam[i]!=-999){
      D->x[i] -= 0.25*(1.0+f[dam[i]]);
    }
    if(sire[i]!=-999){
      D->x[i] -= 0.25*(1.0+f[sire[i]]);
    }

    j=i;
    cnt=0;

    while(j>=0){

      sj=sire[j];
      dj=dam[j];

      if(sj!=-999){
        AN[cnt] = sj;
        li[sj] += 0.5*li[j];
        cnt++;
      }

      if(dj!=-999){ 
        AN[cnt] = dj;
        li[dj] += 0.5*li[j];
        cnt++;
      }

      ai += li[j]*li[j]*D->x[j];

      j=-1;

      for(k=0; k<cnt; k++){   // find eldest individual
        if(AN[k]>j){
          j = AN[k];
        }
      }
      for(k=0; k<cnt; k++){   // delete duplicates
        if(AN[k]==j){
          AN[k] -= (nAP[0]+1); 
        }
      }
    }
    f[i] = ai-1.0;
    for(k=0; k<nAP[0]; k++){
      li[k]  = 0.0;            // reset l to zero except l_ii =1
    }
  }


  for(i=0; i<nAP[0]; i++){  // iterate through each row of l 
      dii[i] = D->x[i];
      D->x[i] = 1.0/(D->x[i]);
  }

  tTD = cs_multiply(tTinv, D);
  tA = cs_multiply(tTD, Tinv);
  A = cs_transpose(tA, TRUE);

  for (i = 0 ; i < A->nzmax; i++){
    iAP[i] = A->i[i];
    xAP[i] = A->x[i];
  }
  for (i = 0 ; i <= A->n; i++){
    pAP[i] = A->p[i];
  }
  nzmaxAP[0] = A->nzmax;

  cs_spfree(Tinv);
  cs_spfree(tTinv);
  cs_spfree(D);
  cs_spfree(tTD);
  cs_spfree(A);
  cs_spfree(tA);
  delete[] AN;
  delete[] li;
}
}
