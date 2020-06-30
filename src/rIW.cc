#include "MCMCglmmcc.h"

extern "C"{  

void rIW(
        double *nu,      // df
        double *V,       // inverse scale matrix
	double *CMP,      // conditional submatrix 
        int *dimP,   // dimension
        int *splitP,   // dimension
        int *n,          // number of samples
        double *output
){         

  GetRNGstate();                                  

  int i, j, dim, split, cnt;
  css *GS;
  cs *G, *Gsamp, *CM;

  dim = dimP[0];
  split = splitP[0];

  G = cs_spalloc(dim, dim, dim*dim, true, false);
	if(split==-999){
	  CM = cs_spalloc(1, 1,1, true, false);
	}else{
      CM = cs_spalloc(dim-split, dim-split, (dim-split)*(dim-split), true, false);
	}
	
  cnt=0;
  for (i = 0 ; i < dim; i++){
    G->p[i] = i*dim;
    for (j = 0 ; j < dim; j++){
      G->i[cnt] = j;
      G->x[cnt] = V[cnt];
      cnt++;
    }
  }
  G->p[dim] = dim*dim;

  GS = cs_schol(0, G);  

  cnt=0;
  if(split==-999){  /* standard inverse wishart */
     for (i = 0 ; i < n[0]; i++){
      Gsamp =  cs_rinvwishart(G, nu[0], GS);
      for (j = 0 ; j < dim*dim; j++){
        output[cnt] = Gsamp->x[j];
        cnt++;
      }
      cs_spfree(Gsamp);
    }
  }else{          /* conditional inverse wishart */
	  cnt=0;
	  for (i = split ; i < dim; i++){
		  CM->p[i-split] = (i-split)*(dim-split);
		  for (j = split ; j < dim; j++){
			  CM->i[cnt] = j-split;
			  CM->x[cnt] = CMP[cnt];
			  cnt++;
		  }
	  }
	  CM->p[(dim-split)] = (dim-split)*(dim-split);
	  cnt=0;
    for (i = 0 ; i < n[0]; i++){
      Gsamp =  cs_rCinvwishart(G, nu[0], split, CM);
      for (j = 0 ; j < dim*dim; j++){
        output[cnt] = Gsamp->x[j];
        cnt++;
      }
      cs_spfree(Gsamp);
    }
  }

  PutRNGstate();
        

  cs_spfree(G);
  cs_spfree(CM);
  cs_sfree(GS);

}
}

