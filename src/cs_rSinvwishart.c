#include "MCMCglmm.h"
#ifndef  Ctol
#define Ctol  1e-12
#endif

cs *cs_rSinvwishart(const cs *A, double nu, int split){

    cs  *A11, *IW, *IW11;
    css *As;
    int nA = A->n,
        nC = nA-split,
        i,
        j,
        cnt = 0;

    A11 = cs_spalloc (split, split, split*split, 1, 0);	
    IW = cs_spalloc (nA, nA, nA*nA, 1, 0);	

    for (i = 0 ; i < split; i++){
      A11->p[i] = i*split;
      for (j = 0 ; j < split; j++){
        A11->i[cnt] = j;
        A11->x[cnt] = A->x[i*nA+j];
        cnt++;
      }
    }

    A11->p[split] = split*split;
	
    As = cs_schol(0, A11);

    IW11 = cs_rinvwishart(A11, nu, As);

    cnt = 0;

    for (i = 0 ; i < split; i++){
      IW->p[i] = i*nA;
      for (j = 0 ; j < split; j++){
        IW->i[cnt] = j;
        IW->x[cnt] = IW11->x[i*split+j];
        cnt++;
      }
      for (j = 0; j < nC; j++){
        IW->i[cnt] = j+split;
        IW->x[cnt] = 0.0;
        cnt++;
      }
    }
    for (i = 0 ; i < nC; i++){
      IW->p[(i+split)] = (i+split)*nA;
      for (j = 0; j < split; j++){
        IW->i[cnt] = j;
        IW->x[cnt] = 0.0;
        cnt++;
      }
      for (j = 0 ; j < nC; j++){
        IW->i[cnt] = j+split;
        if(i==j){
        IW->x[cnt] = 1.0;
        }else{
        IW->x[cnt] = 0.0;
        }
        cnt++;
      }
    } 
    IW->p[nA] = nA*nA;

    cs_spfree(A11);
    cs_spfree(IW11);
    cs_sfree(As);
    return (cs_done (IW, NULL, NULL, 1)) ;	/* success; free workspace, return C */

}


