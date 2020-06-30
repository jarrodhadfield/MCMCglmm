#include "MCMCglmm.h"
# define Dtol  1e-7

cs *cs_rR(const cs *A, double nu, double nuR, const css *As, const cs *Roldinv, double Roldldet, const cs *pG){
    
	cs *Rnew, *Rnewinv, *Ainv;
	double Rnewldet, MH;
        int dimG = A->n;
	int cnt = 0;
	int i, j;
	
	Rnewinv = cs_spalloc (dimG, dimG, dimG*dimG, 1, 0);
	
	for (i = 0 ; i < dimG; i++){
	  Rnewinv->p[i] = i*dimG;
	  for (j = 0 ; j < dimG; j++){
		Rnewinv->i[cnt] = j;
		Rnewinv->x[cnt] = 0.0;
                A->x[i*dimG+j] -= pG->x[i*dimG+j];
 		cnt++;
	  }
	}
	Rnewinv->p[dimG] = dimG*dimG;
		
	cs_cov2cor(A);
	Ainv = cs_inv(A);
	
	Rnew = cs_rinvwishart(Ainv, nu, As);	
	cs_cov2cor(Rnew);
		
	Rnewldet = log(cs_invR(Rnew, Rnewinv));

/*****************************************************/
/*       From Eq A.4 in Liu and Daniels (2006)       */
/*       using \pi_{1} = Eq 6 in Barnard (2000)      */
/*  using \pi_{2} = Eq 3.4 in Liu and Daniels (2006) */
/*****************************************************/

        MH = Roldldet-Rnewldet;
 
	for (i = 0 ; i < dimG; i++){
          MH += log(Roldinv->x[i*dimG+i]);
          MH -= log(Rnewinv->x[i*dimG+i]);
	}

	MH *= 0.5*nuR;

	if(MH<log(runif(0.0,1.0)) || Rnewldet<log(Dtol)){
	  Rnewldet = cs_invR(Roldinv, Rnew);	// save old R	
        }

        for (i = 0 ; i < dimG; i++){
          for (j = 0 ; j < dimG; j++){
 	    Rnew->x[i*dimG+j] *= sqrt((pG->x[i*dimG+i])*(pG->x[j*dimG+j]));
          }
        }

        cs_spfree(Rnewinv);
        cs_spfree(Ainv);

    return (cs_done (Rnew, NULL, NULL, 1)) ;	/* success; free workspace, return C */

}


