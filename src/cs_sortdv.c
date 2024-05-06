#include "MCMCglmm.h"

void cs_sortdv(const cs *A){

    int i;
    double *k= malloc((A->m)*sizeof(double));

    if(k == NULL) {
         Rf_error("could not allocate memory on cs_sortdv");
    }

    for (i = 0 ; i < A->m; i++){   
       k[A->i[i]] = A->x[i]; 
    }
    for (i = 0 ; i < A->m; i++){   
       A->i[i] = i;
       A->x[i] = k[i];
    }

    free(k);
}



