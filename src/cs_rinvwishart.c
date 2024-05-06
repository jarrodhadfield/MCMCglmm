#include "MCMCglmm.h"

cs *cs_rinvwishart(const cs *A, double nu, const css *As){
    
    int m, i, j, cnt;
    cs *T, *IW, *C, *W, *tC;
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
    if(U==NULL){
      PutRNGstate();
      Rf_error("ill-conditioned cross-product: can't form Cholesky factor\n");
    }

    C = cs_multiply(U->L,T);              // t(T)%*%chol(A)
    tC = cs_transpose(C, TRUE);            // t(CI)
    W  = cs_multiply(C,tC);   
    IW = cs_inv(W);                       // crossprod(t(CI))
    cs_spfree(T);
    cs_nfree(U);
    cs_spfree(C);
    cs_spfree(tC);
    cs_spfree(W);

    return (cs_done (IW, NULL, NULL, 1)) ;	/* success; free workspace, return C */

}


