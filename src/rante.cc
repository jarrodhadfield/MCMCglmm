#include "MCMCglmmcc.h"

extern "C"{  

void rante(
        double *locationP,  // effects   
        int *dimP,          // dimension
        int *nlGRP,          // dimension
        int *nkP,            // order of ante
        int *nP,             // number of simulations
        bool *cbetaP,
        bool *cvarP,
        int *iAP,           // inverseA       
	int *pAP,	         
	double *xAP,	    
        int *nzmaxP,      // number of non-zero's inverseA 
        double *GinvP
){         

  int dimG = dimP[0];
  int nlGR = nlGRP[0];
  int nk = nkP[0];
  int n = nP[0];
  bool cbeta = cbetaP[0];
  bool cvar = cvarP[0];
  int Aterm;
  double *ivar = new double[dimG];

  int i, j, cnt, nbeta, start;

  if(nzmaxP[0]==0){
    Aterm = -1;
  }else{
    Aterm = 0;
  }

  cs *Ainv;
  cs *location, *G, *pmuAnte, *pvAnte, *pG;

  location = cs_spalloc(dimG*nlGR, 1, dimG*nlGR, true, false);
  start = 0;

  if(Aterm==0){
    Ainv = cs_spalloc(nlGR, nlGR, nzmaxP[0], true, false);
    for (i = 0; i < nzmaxP[0]; i++){
      Ainv->i[i] = iAP[i];
      Ainv->x[i] = xAP[i];
    }
    for (i = 0; i <= nlGR; i++){
       Ainv->p[i] = pAP[i];
    }
  }

  cnt = 0;
  for (i = 0 ; i < dimG*nlGR ; i++){
    location->i[i] = i;
    location->x[i] = locationP[i];
  }
  location->p[0] = 0; 
  location->p[1] = dimG*nlGR;

  for (i = 0 ; i < dimG; i++){
    ivar[i] = 1.0;
  }


    if(cbeta==1){
      nbeta = nk;
    }else{
      nbeta = dimG*nk-nk*(nk+1)/2;
    }

    pmuAnte = cs_spalloc (nbeta, 1, nbeta, 1, 0);
    pvAnte = cs_spalloc (nbeta, nbeta, nbeta*nbeta, 1, 0);

    cnt = 0;
    for (i = 0 ; i < nbeta; i++){
       pmuAnte->i[i] = i;
       pmuAnte->x[i] = 0.0;
       pvAnte->p[i] = cnt;
       for (j = 0 ; j < nbeta ; j++){
         pvAnte->i[cnt] = j;
         if(i==j){
           pvAnte->x[cnt] = 10e-10;
         }else{
           pvAnte->x[cnt] = 0.0;
         }
         cnt++;
       }
    }
    pmuAnte->p[0]=0;
    pmuAnte->p[nbeta]=nbeta;
    pvAnte->p[nbeta]=nbeta*nbeta;

    pG = cs_spalloc(dimG, dimG, dimG*dimG, true, false);

    cnt = 0;
    for (i = 0 ; i < dimG; i++){
       pG->p[i] = cnt;
       for (j = 0 ; j < dimG; j++){
         pG->i[cnt] = j;
         pG->x[cnt] = 0.0;
         cnt++;
       }
    }
    pG->p[dimG]=dimG*dimG;

    double pnG=0.0;


  GetRNGstate();                                  

  cnt=0;
  for (i = 0 ; i < n; i++){
    G = cs_rAnte(location, start, dimG, nlGR, nk, pmuAnte, pvAnte, Ainv, Aterm, ivar, cvar, pG, pnG);
    for (j = 0 ; j < dimG*dimG ; j++){
      GinvP[cnt] = G->x[j];
      cnt++;
    }
    cs_spfree(G);
  }


  PutRNGstate();
   
  cs_spfree(location);
  cs_spfree(pmuAnte);
  cs_spfree(pvAnte);
  cs_spfree(pG);

  delete [] ivar;

  if(Aterm){
    cs_spfree(Ainv);
  }

}
}

