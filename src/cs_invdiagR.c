#include "MCMCglmm.h"

/* replaces A with the inverse of a diagonalised C matrix and returns the determinant of A*/

double cs_invdiagR(const cs *C, const cs *A){  

    int n, i, j, cnt;
    double det;
	
    n = C->n;
    det=1.0;

    cnt=0;	
    for(i = 0; i<n; i++){  /* diagonalise C before inverting*/
        for(j = 0; j<n; j++){
          if(i==j){
             A->x[cnt] = 1.0/C->x[cnt];
             det  *= A->x[cnt];
          }else{
            A->x[cnt] = 0.0;
          }
          cnt++;  
        }     
    }

    return (det) ;	/* success; free workspace, return C */

}

