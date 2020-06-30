#include "MCMCglmm.h"

cs *cs_dense0(int n, int m){

    int i,j,cnt;
    cs *C;
  
    C = cs_spalloc (n, m, n*m, 1, 0) ;	 /* allocate result */
    if (!C ) return (cs_done (C, NULL, NULL, 0));  

    cnt=0;
    
    for (i = 0 ; i < m; i++){
      C->p[i] = cnt;
      for (j = 0 ; j < n;j++){
        C->i[cnt] = j;
        C->x[cnt] = 0.0;
        cnt++;
      }
    }
    C->p[m] = m*n;
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}



