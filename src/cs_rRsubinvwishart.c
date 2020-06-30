#include "MCMCglmm.h"
#ifndef  Ctol
#define Ctol  1e-12
#endif

cs *cs_rRsubinvwishart(const cs *A, double nu, int split, double nuR, const cs *pG, const cs *oldCM){

    cs  *A22, *Ainv, *IW, *CM, *CMpg, *oldCMinv;
    css *A22s;

    int nA = A->n,
        nC = nA-split,
        i,
        j,
        cnt = 0;
    double Roldldet;

    A22 = cs_spalloc (nC, nC, nC*nC, 1, 0);
    CMpg = cs_spalloc (nC, nC, nC*nC, 1, 0);	
    oldCMinv = cs_spalloc (nC, nC, nC*nC, 1, 0);					

    cnt = 0;
    for (i = 0 ; i < nC; i++){
      A22->p[i] = i*nC;
      CMpg->p[i] = i*nC;
      oldCMinv->p[i] = i*nC;
      for (j = 0 ; j < nC; j++){
        A22->i[cnt] = j;
        A22->x[cnt] = A->x[(i+split)*nA+(j+split)];
        CMpg->i[cnt] = j;
        CMpg->x[cnt] = pG->x[(i+split)*nA+(j+split)];
        oldCMinv->i[cnt] = j;
        oldCMinv->x[cnt] = 1.0;
        cnt++;
      }
    }
    A22->p[nC] = nC*nC;
    CMpg->p[nC] = nC*nC;
    oldCMinv->p[nC] = nC*nC;

    Roldldet = log(cs_invR(oldCM, oldCMinv)); 
    A22s = cs_schol(0, A22);

    CM = cs_rR(A22, nu-split, nuR, A22s, oldCMinv, Roldldet, CMpg);

    Ainv = cs_inv(A);

    IW = cs_rCinvwishart(Ainv, nu, split, CM);

    for (i = 0 ; i < nC*nC; i++){
      oldCM->x[i] = CM->x[i];  // overwite old conditional matrix
    }

    cs_spfree(A22);
    cs_spfree(CM);
    cs_spfree(CMpg);
    cs_spfree(oldCMinv);
    cs_spfree(Ainv);
    cs_sfree(A22s);
    return (cs_done (IW, NULL, NULL, 1)) ;	/* success; free workspace, return C */

}


