#include "MCMCglmmcc.h"

extern "C"{  

/**************************/
/* Meuwissen and Luo 1992 */
/**************************/

void rbv(
        int *id,      
        int *dam,       
        int *sire,         
        double *d,
        double *rbv,     
        int *nidP,
        int *GdimP,        
        double *GinvP,
        int *pedigree,
        int *ggroups,
        double *gmeans,
        int *ngroupP,
        int *nA
){         

  int     i, j, k, cnt, sj, dj;
  int nid = nidP[0];
  int ngroup = ngroupP[0];
  double  *f = new double[nid];
  double  ai;
  double  *AN = new double[nA[0]];
  double  *li = new double[nA[0]];

  int dimG = GdimP[0];

  cs *Ginv, *Grv;
  csn *GinvL;
  css *GinvS;

  Ginv = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
  Grv = cs_spalloc(1, dimG, dimG, true, false);

  cnt = 0;
  for(i = 0 ; i < dimG ; i++){
    Ginv->p[i] = i*dimG;    
    Grv->p[i] = i;
    Grv->i[i] = 0;
    Grv->x[i] = 1.0;
    for(j = 0 ; j < dimG; j++){
      Ginv->i[cnt] = j;
      Ginv->x[cnt] = GinvP[cnt];
      cnt++;
    }
  }
  Ginv->p[dimG] = dimG*dimG;
  Grv->p[dimG] = dimG;

  GinvS = cs_schol(0, Ginv);                 // Symbolic factorisation of G
  GinvL = cs_chol(Ginv, GinvS);              // cholesky factorisation of G^{-1} for forming N(0, G \otimes I)

  GetRNGstate();
  
  if(pedigree[0]==1){
    
    for(i=0; i<nA[0]; i++){
       li[i]=0.0;               // set l to zero
    }
    for(i=0; i<nA[0]; i++){
       AN[i]=-1;               // set AN to zero
    }
         
    for(i=0; i<nid; i++){  // iterate through each row of l 

      li[i] = 1.0;                   // set l_ii to one
      ai=0.0;                        // set a_ii to zero

      if(dam[i]!=-999){
        d[i] -= 0.25*(1.0+f[dam[i]]);
      }
      if(sire[i]!=-999){
        d[i] -= 0.25*(1.0+f[sire[i]]);
      }

      for(j=0; j<dimG; j++){
         Grv->x[j] = rnorm(0.0,sqrt(d[i]));
      }
                  
      cs_ltsolve(GinvL->L, Grv->x);

      for(j=0; j<dimG; j++){
        rbv[j*nid+i] = Grv->x[j];     
      }

      if(sire[i]!=-999){  
        for(j=0; j<dimG; j++){
          rbv[j*nid+i] += 0.5*(rbv[j*nid+sire[i]]);
        }  
      }else{
        for(j=0; j<dimG; j++){
          rbv[j*nid+i] += 0.5*gmeans[ggroups[i]+j*ngroup];
        } 
      }
      if(dam[i]!=-999){  
        for(j=0; j<dimG; j++){
           rbv[j*nid+i] += 0.5*(rbv[j*nid+dam[i]]);
        }  
      }else{
        for(j=0; j<dimG; j++){
          rbv[j*nid+i] += 0.5*gmeans[ggroups[i]+j*ngroup];
        } 
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

        ai += li[j]*li[j]*d[j];

        j=-1;
        for(k=0; k<cnt; k++){   // find eldest individual
          if(AN[k]>j){
            j = AN[k];
          }
        }
        for(k=0; k<cnt; k++){   // delete duplicates
          if(AN[k]==j){
            AN[k] -= (nid+1); 
          }
        }
      }
      f[i] = ai-1.0;
      for(k=0; k<nid; k++){
        li[k]  = 0.0;            // reset l to zero except l_ii =1
      }
    }

  }else{

    for(i=0; i<nid; i++){  // iterate through each row of l 

      for(j=0; j<dimG; j++){
        Grv->x[j] = rnorm(0.0,sqrt(d[i]));
      }
                  
      cs_ltsolve(GinvL->L, Grv->x);

      for(j=0; j<dimG; j++){
        rbv[j*nid+i] = Grv->x[j];     
      }
      if(dam[i]!=-999){  
        for(j=0; j<dimG; j++){
          rbv[j*nid+i] += rbv[j*nid+dam[i]];
        }  
      }
    }  
  }

  PutRNGstate();

  cs_spfree(Ginv);
  cs_spfree(Grv);
  cs_sfree(GinvS);                
  cs_nfree(GinvL);        
  delete[] f;
  delete[] AN;
  delete[] li;
}
}
