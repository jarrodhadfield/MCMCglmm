#include "MCMCglmm.h"

cs *cs_rwishart(const cs *A, double nu, const css *As){
    
    int m, i, j, cnt;
    cs *T, *C, *W, *tC;
    csn *U;
    m = A->n;

    T = cs_spalloc (m, m, m*(m+1)/2, 1, 0) ;	 
    if (!T ) return (cs_done (T, NULL, NULL, 0));   

    double df = nu;
    cnt = 0;

    for(i = 0; i<m; i++){
      T->p[i] = cnt;  
      T->i[cnt] = i;    
      T->x[cnt] = sqrt(rchisq(df));
      cnt++;
      for(j = i+1; j<m; j++){
        T->i[cnt] = j;
        T->x[cnt] = rnorm(0.0,1.0);
        cnt++;
      } 
      df--;
    }
    T->p[m] = m*(m+1)/2;

    U = cs_chol(A, As);  
    C = cs_multiply(U->L,T);              // t(T)%*%chol(A)
    tC = cs_transpose(C, TRUE);            // t(CI)
    W  = cs_multiply(C,tC);   
    cs_spfree(T);
    cs_nfree(U);
    cs_spfree(C);
    cs_spfree(tC);

    return (cs_done (W, NULL, NULL, 1)) ;	/* success; free workspace, return C */

}


