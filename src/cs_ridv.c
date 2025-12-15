#include "MCMCglmm.h"


cs *cs_ridv(const cs *A, double nu, int split, const cs *CM){
    
	cs *Vnew;
    int nA = A->n;
    int nC = nA-split;
	int cnt = 0;
	int i, j;
	double ss = 0.0;

	Vnew = cs_spalloc (nA, nA, nA*nA, 1, 0);
	
    if(split<0){

		for (i = 0 ; i < nA; i++){
	        ss += A->x[i*nA+i];
	    }    

	    ss = ss/rchisq(nA*nu);

		for (i = 0 ; i < nA; i++){
			Vnew->p[i] = i*nA;
			for (j = 0 ; j < nA; j++){
				Vnew->i[cnt] = j;
				if(i==j){
				  Vnew->x[cnt] = ss;
				}else{	
				  Vnew->x[cnt] = 0.0;
				}  
				cnt++;
			}
		}
	}else{
		
		for (i = 0 ; i < split; i++){
			ss += A->x[i*nA+i];
		}    

		ss = ss/rchisq(split*nu);

		for (i = 0 ; i < split; i++){
			Vnew->p[i] = i*nA;
			for (j = 0 ; j < split; j++){
				Vnew->i[cnt] = j;
				if(i==j){
				Vnew->x[cnt] = ss;
				}else{	
				Vnew->x[cnt] = 0.0;
				}  
				cnt++;
			}
			for (j = 0; j < nC; j++){
				Vnew->i[cnt] = j+split;
				Vnew->x[cnt] = 0.0;
				cnt++;
			}
		}
		for (i = 0 ; i < nC; i++){
			Vnew->p[(i+split)] = (i+split)*nA;
			for (j = 0; j < split; j++){
				Vnew->i[cnt] = j;
				Vnew->x[cnt] = 0.0;
				cnt++;
			}
			for (j = 0 ; j < nC; j++){
				Vnew->i[cnt] = j+split;
				Vnew->x[cnt] = CM->x[i*nC+j];
				cnt++;
			}
		} 
	}	
    Vnew->p[nA] = nA*nA;
		
    return(cs_done (Vnew, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}

