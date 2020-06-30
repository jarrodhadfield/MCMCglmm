#include "MCMCglmm.h"

cs *cs_rAnte(const cs *location, int start, int dimG, int nlGR, int nk, const cs *pmuAnte, const cs *pvAnte, const cs *Ainv, int Aterm, double *ivar, int cvar, const cs *pG, double pnG){

    int i, j, k, cnt, cnt2, triangle, nbeta, cbeta;
    cs *X, *tX, *cholGinv, *tcholGinv, *Ginv, *Rinv, *beta_star, *z_star, *tXRinv, *tXRinvX, *MME, *beta, *beta_tmp, *pred, *G, *location_tmp, *tlAl, *tlA, *tl; 
    csn *pvAnteL, *L, *RinvL;
    css *pvAnteS, *S, *RinvS;

    triangle = nk*(nk+1)/2;
    nbeta = pmuAnte->m;

    if(nbeta==nk){
      cbeta=1;
    }else{
      cbeta=0;
    }

    pvAnteS = cs_schol(0, pvAnte);  // really, this just needs to be done once
    pvAnteL = cs_chol(pvAnte, pvAnteS); 

    if(cbeta==1){
       X = cs_spalloc (nlGR*dimG, nk, nlGR*(dimG*nk-triangle), 1, 0);
       cnt = 0;
       for(k=1; k<=nk; k++){ 
         X->p[k-1] = cnt; 
         for(j=k; j<dimG; j++){ 
           for(i=0; i<nlGR; i++){
              X->i[cnt] = j*nlGR+i;    
              X->x[cnt] = location->x[start+(j-k)*nlGR+i];
              cnt++;
           }
         }
       }
       X->p[nbeta] = nlGR*(dimG*nk-triangle); 
    }else{
       X = cs_spalloc (nlGR*dimG, nbeta, nlGR*nbeta, 1, 0) ;
       cnt = 0;
       cnt2 = 0;
       for(k=1; k<=nk; k++){ 
         for(j=k; j<dimG; j++){ 
           X->p[cnt2] = cnt; 
           cnt2++;
           for(i=0; i<nlGR; i++){
              X->i[cnt] = j*nlGR+i;    
              X->x[cnt] = location->x[start+(j-k)*nlGR+i];
              cnt++;
           }
         }
       }
       X->p[nbeta] = nlGR*nbeta; 
    }	 
    tX = cs_transpose(X, TRUE);                               

    // sample regression coefficients from their prior

    beta_star = cs_spalloc (nbeta, 1, nbeta, 1, 0);
    beta = cs_spalloc (nbeta, 1, nbeta, 1, 0);

    beta_star->p[0]=0;
    beta_star->p[1]=nbeta;
    beta->p[0]=0;
    beta->p[1]=nbeta;

    for(i=0; i<nbeta; i++){
      beta_star->i[i] = i;
      beta_star->x[i] = rnorm(0.0,1.0);
      beta->i[i] = i;
      beta->x[i] = 0.0;
    }                   

    cs_ltsolve(pvAnteL->L, beta_star->x);

    for(i=0; i<nbeta; i++){
       beta_star->x[i]  += pmuAnte->x[i];
    }                        
                       
    // create solve(R) and sample e_star-y

    location_tmp = cs_spalloc (nlGR, dimG, nlGR*dimG, 1, 0);

    cnt=0;
    for (i = 0 ; i < dimG; i++){
      location_tmp->p[i] = cnt;
      for (j = 0 ; j < nlGR; j++){
        location_tmp->i[cnt] = j;
        cnt++;
      }
    }
    location_tmp->p[dimG] = nlGR*dimG;


    z_star = cs_spalloc (nlGR*dimG, 1, nlGR*dimG, 1, 0);

    
    for (i = 0 ; i < dimG; i++){  
        ivar[i] = 1.0/ivar[i];
    }

    if(Aterm>=0){
      Rinv = cs_kroneckerDA(ivar,dimG, Ainv);
      RinvS = cs_schol(1, Rinv);
      RinvL = cs_chol(Rinv, RinvS);

      if(RinvL==NULL){
        error("problems with ginverse in antependence model");
      }
      
// Jarrod: could speed things up by using AinvL and AinvS and going through one at a time

      cnt = 0;
      for (i = 0 ; i < dimG; i++){
        for (j = 0 ; j < nlGR; j++){
          z_star->i[cnt] = cnt;
          location_tmp->x[cnt] = rnorm(0.0,1.0);  
          cnt++;
        }
      }
      z_star->p[0] = 0;
      z_star->p[1] = dimG*nlGR;
        
      cs_ltsolve(RinvL->L,  location_tmp->x);
      for(i=0; i<(nlGR*dimG); i++){  
        z_star->x[i]  = location_tmp->x[RinvS->pinv[i]]-location->x[start+i];
      }
        
      cs_nfree(RinvL);
      cs_sfree(RinvS);
        
    }else{
      Rinv = cs_kroneckerDI(ivar,dimG, nlGR);
      cnt = 0;
      for (i = 0 ; i < dimG; i++){
        for (j = 0 ; j < nlGR; j++){
          z_star->i[cnt] = cnt;
          z_star->x[cnt] = rnorm(0.0,sqrt(1.0/ivar[i]))-location->x[start+cnt];
          cnt++;
        }
      }
      z_star->p[0] = 0;
      z_star->p[1] = dimG*nlGR;
    }


    // create t(X)%*%solve(R)%*%X+solve(B)

    tXRinv = cs_multiply(tX, Rinv);
    tXRinvX = cs_multiply(tXRinv, X);
    MME = cs_add(tXRinvX, pvAnte, 1.0, 1.0); 

    // create X%*%beta_star+e_star-y

    cs_gaxpy(X, beta_star->x, z_star->x); 

    // create t(X)%*%solve(R)%*%(X%*%beta_star+e_star-y)

    beta_tmp = cs_multiply(tXRinv, z_star);  

    for (i = 0 ; i < nbeta; i++){
      beta->x[beta_tmp->i[i]] = -beta_tmp->x[i];
    }
    for (i = 0 ; i < nbeta; i++){
      beta_tmp->x[i] = 0.0;
    }

    S = cs_schol(1, MME);
    L = cs_chol(MME, S); 

    if(L==NULL){
      PutRNGstate();
      error("antedependence equations singular: use a (stronger) prior\n");
    }

    cs_ipvec (S->pinv, beta->x, beta_tmp->x, MME->n);	 // x = P*b 
    cs_lsolve(L->L, beta_tmp->x);                        // x = L\x 
    cs_ltsolve (L->L, beta_tmp->x);		         // x = L'\x 
    cs_pvec (S->pinv, beta_tmp->x, beta->x, MME->n);     // b = P'*x 

    for (i = 0 ; i < nbeta ; i++){
      beta->x[i] += beta_star->x[i];
    }

    pred = cs_multiply(X, beta);

    for(i=0; i<(nlGR*dimG); i++){
      location_tmp->x[i] = location->x[start+i];
    }
    for(i=0; i< pred->p[1]; i++){
      location_tmp->x[pred->i[i]] -= pred->x[i];
    }

    tl = cs_transpose(location_tmp, TRUE);

    if(Aterm>=0){
      tlA = cs_multiply(tl, Ainv);
      tlAl = cs_multiply(tlA, location_tmp);
      cs_spfree(tlA);
    }else{
      tlAl = cs_multiply(tl, location_tmp);
    }

    if(cvar==1){
       ivar[0] = 0.0;
       for(i=0; i<dimG; i++){
           ivar[0]  += tlAl->x[i*dimG+i];
       }
       ivar[0] += pG->x[0];
       ivar[0] /= rchisq((double)nlGR*dimG+pnG);
       for (i = 1 ; i < dimG; i++){
         ivar[i] = ivar[0];
       }
    }else{
       for (i = 0 ; i < dimG; i++){
         ivar[i] = (tlAl->x[i*dimG+i]+pG->x[i*dimG+i])/rchisq((double)nlGR+pnG);
       }
    }

    cholGinv = cs_spalloc (dimG, dimG, dimG*dimG, 1, 0);

    cnt=0;

    for (i = 0 ; i < dimG; i++){
       cholGinv->p[i] = cnt;
       for (j = 0 ; j < dimG ; j++){
         cholGinv->i[cnt] = j;
         cholGinv->x[cnt] = 0.0;
         cnt++;
       }
    }
    cholGinv->p[dimG] = dimG*dimG;

    for (i = 0 ; i < dimG; i++){
       cholGinv->x[i*dimG+i] = 1.0/sqrt(ivar[i]);
    }

    if(cbeta==1){
      for(j=1; j <= nk; j++){
        for(i=0; i<(dimG-j); i++){ 
          cholGinv->x[i*dimG+i+j] = -beta->x[j-1]/sqrt(ivar[i+j]);
        }
      }
    }else{
      cnt=0;
      for(j=1; j <= nk; j++){
        for(i=0; i<(dimG-j); i++){ 
          cholGinv->x[i*dimG+i+j] = -beta->x[cnt]/sqrt(ivar[i+j]);
          cnt++;
        }
      }
    }

    tcholGinv = cs_transpose(cholGinv, TRUE);
    Ginv = cs_multiply(tcholGinv, cholGinv);

    G = cs_inv(Ginv);

    cs_spfree(X);
    cs_spfree(tX);
    cs_spfree(Rinv);
    cs_spfree(beta_star);
    cs_spfree(z_star);
    cs_spfree(tXRinv);
    cs_spfree(tXRinvX);
    cs_spfree(MME);
    cs_spfree(beta);
    cs_spfree(beta_tmp);
    cs_spfree(pred); 
    cs_spfree(cholGinv);
    cs_spfree(tcholGinv);
    cs_spfree(location_tmp);
    cs_spfree(Ginv);
    cs_spfree(tl);
    cs_spfree(tlAl); 
    cs_nfree(pvAnteL);
    cs_nfree(L);
    cs_sfree(pvAnteS);
    cs_sfree(S);
    
    return (cs_done (G, NULL, NULL, 1)) ;	/* success; free workspace, return C */

}


