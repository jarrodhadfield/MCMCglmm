#include "MCMCglmm.h"
#ifndef  Ctol
#define Ctol  1e-12
#endif

cs *cs_rCinvwishart(const cs *A, double nu, int split, const cs *CM){

    cs  *T1inv, *T1, *A11, *A12, *A21, *A22, *A11Schur, *A11Schurinv, *halfA11Schur, *CinvSchur, *A11inv, *varT2inv, *IW, *IW11, *halfA11SchurT, *halfCSchur;
    csn *varT2invL;
    css *As, *T2s;
    int nA = A->n,
        nC = nA-split,
        i,
        j,
        cnt = 0,
        cnt2 = 0;
    double Rv[nC*split];

    A11 = cs_spalloc (split, split, split*split, 1, 0);	
    A12 = cs_spalloc (nC, split, nC*split, 1, 0);	
    A22 = cs_spalloc (nC, nC, nC*nC, 1, 0);		
    IW = cs_spalloc (nA, nA, nA*nA, 1, 0);	

    for (i = 0 ; i < split; i++){
      A11->p[i] = i*split;
      A12->p[i] = i*nC;
      for (j = 0 ; j < split; j++){
        A11->i[cnt] = j;
        A11->x[cnt] = A->x[i*nA+j];
        cnt++;
      }
      for (j = 0 ; j < nC; j++){
        A12->i[cnt2] = j;
        A12->x[cnt2] = A->x[i*nA+j+split];
        cnt2++;
      }
    }

    A11->p[split] = split*split;
    A12->p[split] = nC*split;
    A21 = cs_transpose(A12, TRUE);

    A11inv = cs_inv(A11);
   
    cnt = 0;
    for (i = 0 ; i < nC; i++){
      A22->p[i] = i*nC;
      for (j = 0 ; j < nC; j++){
        A22->i[cnt] = j;
        A22->x[cnt] = A->x[(i+split)*nA+(j+split)];
        cnt++;
      }
    }
    A22->p[nC] = nC*nC;
	
    halfA11Schur = cs_multiply(A11inv, A21);
    A11Schur = cs_multiply(A12, halfA11Schur);
	
    for (i = 0 ; i < (nC*nC); i++){
       A11Schur->x[i] =  A22->x[i] - A11Schur->x[i];
    }
	
    As = cs_schol(0, A11);

    T1inv = cs_rinvwishart(A11, nu, As);

    for(i=0; i<(nC*split); i++){
       Rv[i] = rnorm(0.0,1.0);
    }

    T1 = cs_inv(T1inv);
    A11Schurinv = cs_inv(A11Schur);
    varT2inv = cs_kroneckerA(T1, A11Schurinv);


    T2s = cs_schol(0, varT2inv);
    varT2invL = cs_chol(varT2inv, T2s);
    cs_ltsolve(varT2invL->L, Rv);
    
    for (i = 0 ; i < (nC*split); i++){
      halfA11Schur->x[i] +=  Rv[i];  
      halfA11Schur->x[i] *= -1.0;  // This equals -T2 from Koorsgard et al 
    }

    halfA11SchurT = cs_transpose(halfA11Schur, TRUE);  // -T2'
    halfCSchur = cs_multiply(CM, halfA11SchurT);       // -T3^{-1}T2'   
    CinvSchur = cs_multiply(halfA11Schur, halfCSchur); // T2T3^{-1}T2' 

    IW11 = cs_add(T1inv, CinvSchur, 1.0,1.0);          // T1^{-1}+T2T3^{-1}T2'
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
        IW->x[cnt] = halfCSchur->x[i*nC+j];
        cnt++;
      }
    }
    for (i = 0 ; i < nC; i++){
      IW->p[(i+split)] = (i+split)*nA;
      for (j = 0; j < split; j++){
        IW->i[cnt] = j;
        IW->x[cnt] = halfCSchur->x[j*nC+i];
        cnt++;
      }
      for (j = 0 ; j < nC; j++){
        IW->i[cnt] = j+split;
        IW->x[cnt] = CM->x[i*nC+j];
        cnt++;
      }
    } 
    IW->p[nA] = nA*nA;

    cs_spfree(T1inv);
    cs_spfree(A11);
    cs_spfree(A12);
    cs_spfree(A21);
    cs_spfree(A22);
    cs_spfree(A11inv);
    cs_spfree(A11Schur);
    cs_spfree(halfA11Schur);
    cs_spfree(CinvSchur);
    cs_spfree(halfCSchur);
    cs_spfree(halfA11SchurT);
    cs_spfree(T1);
    cs_spfree(A11Schurinv);
    cs_spfree(varT2inv);
    cs_spfree(IW11);
    cs_nfree(varT2invL);
    cs_sfree(T2s);
    cs_sfree(As);
    return (cs_done (IW, NULL, NULL, 1)) ;	/* success; free workspace, return C */

}


