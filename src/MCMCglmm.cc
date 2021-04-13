#include "MCMCglmmcc.h"

extern "C"{  

/***************************************************************************************************/
/* Matrices are stored in compressed-column format where i indicates the row indicies of non-zero  */
/* values, p indexes the first i of each column, and x the actual non-zero values.  For example,if */
/*                                                                                                 */
/*       0 1 0            p = {0, 1, 3, 4}                                                         */
/*   X = 1 0 2    then,   i = {1, 0, 2, 1}                                                         */
/*       0 1 0            x = {1.0, 1.0, 1.0, 2.0}                                                 */
/*                                                                                                 */
/* where the final element of p is always length(i).                                               */
/* dim is a vector with the number of rows and columns, and nzmax the max number of non-zero values*/
/***************************************************************************************************/

void MCMCglmm(
        double *yP,      // read in data
        double *y2P,     // additional data variables (censoring points, truncation points, binomial n etc)
        double *liabP,   // starting liabilities (= y if gaussian, and no missing data)
        int *mvtype,     // missing data types 
        int *nyP,        // number of records
      	int *dimP,	 // dimensions of X,Z and L 
        int *nzmaxP,     // number of non-zero's in X,Z,L and the A's 
        int *iXP,        // X       
      	int *pXP,	         
      	double *xXP,	               
        int *iZP,        // Z       
      	int *pZP,	         
      	double *xZP,	                  
        int *iAP,        // A       
      	int *pAP,	         
      	double *xAP,	               
        int *AtermP,         // integer: each G is associated with kth A (-1 means A=I)
        int *GRdimP,         // GR structure dimensions
        int *levelsRP,       // number of GR levels 
        int *updateP,        // type of VCV update then the dimension G in a G-R VCV (0 = R and G are separable) and the update for G-R VCV.
        int *splitP, 
        int *nGP,            // number of G structures and number of R structures  
        double *GRinvP,      // GR starting values 
        double *GRvpP,       // GR prior V 
        double *GRnpP,       // GR prior df:
        int *chainP,         // Number of MCMC iterations, thinning interval and burnin. Should posterior distribution of random effects and liabilities be stored	
        double *LocP,        // Posterior distribution of MME solutions (fixed and random)
        double *VarP,        // Posterior distribution of (co)variances
        double *PLiabP,      // Posterior distribution of liabilities
        int *familyP,        // distribution of response variables
        double *propCP,      // proposal distribition for liabilities
        bool *verboseP,      // print iterations and MH acceptance ratio to screen
        double *BvpP,        // inverse prior covariance matrix for beta  with final number indicating diagonal (1) or non diagonal (0)      
        double *BmupP,       // prior mean vector for beta     
        int *mfacP,          // vector of J-1 levels for each multinomial response 
      	int *observedP,	     // vector of 1 (observed) and 0 (missing)
      	int *diagP,          // is a us matrix in fact diagonal?
        int *AMtuneP,        // should adaptive Metroplis algorithm be used
      	int *DICP,	     // should DIC be computed
        double *dbarP,
        int *proposal,
        int *ncutpointsP,    // number of cutpoints with -Inf, 0, C1,.... Cncutpoints, Inf
        int *nordinalP,      // number of ordinal traits  
        double *stcutpointsP,// starting vector of cutpoints    
        double *CPP,
        double *AmupP,
        int *iAVpP,          // priors for redundant parameters in parameter expanded models.      
        int *pAVpP,        
        double *xAVpP,     
        int *nzmaxAVpP,        
        int *PXtermsP,
        int *iLXP,              // Design matrix to form Gianola & Sorensen's Lambda         
        int *pLXP,        
        double *xLXP,
        double *lambdaP,        
        double *LvpP,          // prior for structural parameters
        double *LmupP,
        double *nanteBP,       // orders of antedependence beta's per structure then numbers of betas and then whether constant variance
        double *anteBvpP,      // inverse prior covariance matrix for antedpendence betas     
        double *anteBmupP,    // prior mean vector for antedpendence betas
        bool *truncP,          // should latent variables be truncated tp prevent under/overflow for cat/mult/thresh/ord distributions and zero-processes in zi/za/hu?
        int *me_rtermP,        // vector indexing which me_outcome a set of resdiuals are informative about 
        double *me_prior_probP,   //  prior category probabilities for each me_outcome    
        double *me_XiP,       // me interaction design matrix
        int *meP       // number of po terms, where first covariate for each term appears in X, number of covariates in each term

){         

int     i, j, k,l,p,cnt,cnt2,cnt3,cnt4, rterm,itt,record,dimG,nthordinal,

        nG = nGP[0],          // number of G structures
        nR = nGP[1],          // number of R structures 
        nGR = nG+nR,          // number of variance structures
        ny = nyP[0],          // number of records
        nrterms,              // number of r-blocks
        *nlGR = levelsRP,     // number of random levels 
        *GRdim = GRdimP,      // dimensions of variance structures
        ncond=0,
        nkeep=0,
        nitt = chainP[0], 
        thin = chainP[1], 
        burnin = chainP[2], 
        post_cnt = 0,
        tvc =0,          // total number of (co)variance components 
        nordinal = nordinalP[0],
        diagR = diagP[0];  // is the us R-structure diagionalised 0-no, 1-yes, 2-yes with identical variances

// Multinomial counters 
 
int     nthmnl = 0;  // counter for the lth level of nth multinomial

double   mndenom1 = 1.0,  // sum(exp(l_{j}) for all j) and l_old
         mndenom2 = 1.0,  // sum(exp(l_{j}) for all j) and l_new
         u;               // slice sampling runifs

int     nrowX =  dimP[0], ncolX =  dimP[1],  nzmaxX = nzmaxP[0]; 
int     nrowZ =  dimP[2], ncolZ =  dimP[3],  nzmaxZ = nzmaxP[1];
int     nrowLX =  dimP[4], ncolLX =  dimP[5],  nzmaxLX = nzmaxP[2];

int     nL;
bool    path;

if(nrowLX!=0){
  if(nrowLX!=nrowX){
    path = TRUE;         // path model defined within residual blocks
    nL = ncolLX/GRdimP[nG];
  }else{
    path = FALSE;        // path models cross residual blocks
    nL = ncolLX/nrowLX;  // number of structural parameters
  }
}else{
  path = FALSE;
  nL = 0;
}


int 	dimAS =  ncolX+ncolZ;
int     *nMH = new int[nGR];
int     nalpha = 0; 
 
bool    pr = chainP[3];
bool    pl = chainP[4];
bool    cp = nordinal>0;
bool    missing = FALSE;
int    covu = updateP[nGR];


double detLambda[2];   // holds old and proposed fabs(det(Lambda))
int sign_detLambda[2];    // holds sign of fabs(det(Lambda))
int lambda_old = 0;    // indexes freed Lambda
int lambda_new = 1;    // indexes freed Lambda
int *diagLambdaL = new int[nrowLX]; // indexes diagonal elements of L and U from LU decomposition of Lambda
int *diagLambdaU = new int[nrowLX];

double log_alphaL;     // MH ratio for lambda

cnt=0;
nrterms=0;
for(k=nG; k<nGR; k++){
  nMH[k] = 0;
  if(GRdim[k]>cnt){cnt = GRdim[k];}
  for(i=0; i<nlGR[k]; i++){
    if(mvtype[i+nrterms]!=2){
      missing=TRUE;  // indicates whether there are any non-Gaussian data or missing observations
      if(mvtype[i+nrterms]<0){
        nMH[k]++; 
      }
    }
  }
  nrterms+=nlGR[k];
}

int  *present = new int[nrterms+nR]; // vector of 0 and 1's for ztmultinomial

int *cond = new int[cnt];
int *keep = new int[cnt];

for(k=0; k<nGR; k++){           
   nalpha += PXtermsP[k]*GRdim[k];
}

double  densityl1,
        densityl2,
        dbar = 0.0,
        mdbar = 0.0,
        *Eaccl = new double[nGR],
        *alpha_star = new double[nGR],
        rACCEPT = 0.9,
        qACCEPT = 2.0,
        *t = new double[nGR*2],
        *sd = new double[nGR*2],
        *wn = new double[nGR*2],
        *zn = new double[nGR*2],
        *ldet = new double[nR+nG],
        ldet_rr,
        interval,
        remainder,
        logitt = 1e+35, // truncation limits to stop under/overflow 
        probitt = 1e+35;
 	
if(truncP[0]){
  logitt = 20.0; // truncation limits to stop under/overflow 
  probitt = 7.0;
}
//double inf = std::numeric_limits<double>::max();

int *cumsum_ncutpoints = new int[nordinal+1];
    cumsum_ncutpoints[0] = 0;
int ncutpoints_store = 0;
int *cutpoints_updated = new int[nordinal];

if(cp){
  for(i=0; i<nordinal; i++){          
    cumsum_ncutpoints[i+1] = ncutpointsP[i]+cumsum_ncutpoints[i];  
    if(ncutpointsP[i]>3){
      ncutpoints_store += ncutpointsP[i]-3;
    }
  }
} 

double  *oldcutpoints = new double[cumsum_ncutpoints[nordinal]],
        *newcutpoints = new double[cumsum_ncutpoints[nordinal]],
        *sdcp = new double[nordinal],
        *wncp = new double[nordinal],
        *zncp = new double[nordinal],
        *accp = new double[nordinal],
        *cutpointMHR = new double[nordinal];

if(cp){
  cnt = 0;
  for(i=0; i<nordinal; i++){
    sdcp[i] =1.0;
    wncp[i] =1.0;
    zncp[i] =1.0;
    accp[i] =0.0;
    if(truncP[0]){
      stcutpointsP[cnt]=-probitt;
      stcutpointsP[cnt+ncutpointsP[i]-1]=probitt;
    }
    for(j=0; j<ncutpointsP[i]; j++){
      oldcutpoints[cnt] = stcutpointsP[cnt];
      newcutpoints[cnt] = stcutpointsP[cnt];
      cnt ++;
    }
  }
}

if(ncutpoints_store==0){cp=FALSE;}  // even if cutpoints are present they do not need to be updated if there are only 2 categories

for(k=0; k<nGR; k++){
  t[k]=0.0;
  sd[k]=1.0;
  wn[k]=0.0;
  zn[k]=0.0;
  t[k+nGR]=0.0;
  sd[k+nGR]=1.0;
  wn[k+nGR]=0.0;
  zn[k+nGR]=0.0;
  Eaccl[k] = 0.0;
  alpha_star[k] = (-double(GRdimP[k])/(1.0-2.75*double(GRdimP[k])))-0.133;  // optimal acceptance ratio
}

int nme = meP[0];
int me_rterm;
int *stme = new int[nme];
int *nmec = new int[nme];
int *nmei = new int[nme];
int *nmeo = new int[nme];
int *me_outcome = new int[nrowX*nme];
double tprob;

cnt=0;
for(i=0; i<nme; i++){
  stme[i] = meP[1+i];
  nmec[i] = meP[nme+1+i];
  nmei[i] = meP[2*nme+1+i];
  nmeo[i] = meP[3*nme+1+i];
}  

cs      *X, *Z, *W,  *Wt, *KRinv, *WtmKRinv, *WtmKRinvtmp, *M, *Omega, *MME, *zstar, *astar, *astar_tmp, *location, *location_tmp, *linky, *mulinky,  *pred, *mupred, *dev, *bvA, *bvAbv, *tbv, *pvB, *pmuB, *Brv, *Xalpha, *tXalpha, *Alphainv, *muAlpha, *XtmKRinv, *XtmKRinvtmp, *alphaM, *alphaMME, *alphaastar, *alphapred, *alphazstar, *alphaastar_tmp, *alphalocation, *alphalocation_tmp, *Worig, *LambdaX, *pvL, *pmuL, *Lrv, *I, *linky_orig, *Y, *tY, *ILY, *w, *tYKrinv, *tYKrinvY, *tYKrinvw, *lambda_dev, *tl, *tlV, *tlVl, *exLambda, *exLambdaX, *mulambda, *G_rr, *Ginv_rr, *beta_rr;

csn	*L, *pvBL, *alphaL, *AlphainvL, *pvLL, *tYKrinvYL;
css     *S, *pvBS, *alphaS, *AlphainvS, *pvLS, *LambdaS, *tYKrinvYS, *GinvS_rr;

cs*     *Ginv = new cs*[nGR];
cs*     *muG = new cs*[nGR];	
cs*     *propC = new cs*[nGR*2];
cs*     *propCinv = new cs*[nGR*2];
cs*     *muC = new cs*[nGR*2];
cs*     *linki = new cs*[nGR];
cs*     *linki_tmp = new cs*[nGR];
cs*     *linki_tmp2 = new cs*[nGR];
cs*     *predi = new cs*[nGR];
cs*     *G = new cs*[nGR];
cs*     *pG = new cs*[nGR];
cs*     *CM = new cs*[nGR];	
cs*     *Gtmp = new cs*[nGR];
cs*     *Grv = new cs*[nGR];
css*    *GinvS = new css*[nGR];
csn*    *GinvL = new csn*[nGR];
css*    *propCinvS = new css*[nGR*2];
csn*    *propCinvL = new csn*[nGR*2];
cs*     *KGinv = new cs*[nGR];
csn*    *KGinvL = new csn*[nGR];
css*    *KGinvS = new css*[nGR];
cs*     *lambda = new cs*[2];
cs*     *lambdaI = new cs*[2];
cs*     *Lambda = new cs*[2];
csn*    *LambdaLU = new csn*[2];
cs*     *Lambda_tmp = new cs*[2];
cs*     *A = new cs*[nGR];
cs*     *bv = new cs*[nGR];
cs*     *bv_tmp = new cs*[nGR];
cs*     *pvAnte = new cs*[nGR];
cs*     *pmuAnte = new cs*[nGR];
cs*     *ivar = new cs*[nGR];
cs*     *me_prior_prob = new cs*[nme];
cs*     *me_post_prob = new cs*[nme];
cs*     *me_Xi = new cs*[nme];

//Rprintf("read X");
/****************************************************/ 
/* read in fixed-effects design matrix X and priors */
/****************************************************/

// read in fixed-effects design matrix X 

X = cs_spalloc(nrowX, ncolX, nzmaxX, true, false); /* X =  0  X_2  0  */
                                                   /*      0   0  X_3 */
for (i = 0 ; i < nzmaxX ; i++){
  X->i[i] = iXP[i];
  X->x[i] = xXP[i];
}
for (i = 0 ; i <= ncolX ; i++){
  X->p[i] = pXP[i];
}

// read in prior covaraince matrix and mean vector for fixed effects.

pvB = cs_spalloc(ncolX, ncolX, pow(ncolX, 2.0-BvpP[ncolX*ncolX]), true, false);
pmuB = cs_spalloc(ncolX, 1, ncolX, true, false);
Brv = cs_spalloc(ncolX, 1, ncolX, true, false);

cnt = 0;
for (i = 0 ; i < ncolX ; i++){
 pmuB->i[i] = i;
 Brv->i[i] = i;
 pmuB->x[i] = BmupP[i];
}
pmuB->p[0] = 0; 
pmuB->p[1] = ncolX;
Brv->p[0] = 0; 
Brv->p[1] = ncolX;

if(BvpP[ncolX*ncolX]>0.5){   // prior for beta is diagonal
 for (i = 0 ; i < ncolX ; i++){
   pvB->p[i] = i; 
   pvB->i[i] = i;
   pvB->x[i] = BvpP[ncolX*i+i];
 }
 pvB->p[ncolX] = ncolX;
}else{
 cnt = 0;
 for (i = 0 ; i < ncolX ; i++){
   pvB->p[i] = i*ncolX; 
   for (j = 0 ; j < ncolX ; j++){
     pvB->i[cnt] = j;
     pvB->x[cnt] = BvpP[cnt];
     cnt++;
   }
 }
 pvB->p[ncolX] = ncolX*ncolX;
}

pvBS = cs_schol(0, pvB);                    // Symbolic factorisation of B
pvBL = cs_chol(pvB, pvBS);                  // cholesky factorisation of B^{-1} for forming N(0, B)

/*************************************/ 
/* read in path analytic structures  */
/*************************************/

if(nL>0){

// read in path analytic design matrix 
                                 
  LambdaX = cs_spalloc(nrowLX, ncolLX, nzmaxLX, true, false); 
                                                             
  for (i = 0 ; i < nzmaxLX ; i++){
    LambdaX->i[i] = iLXP[i];
    LambdaX->x[i] = xLXP[i];
  }
  for (i = 0 ; i <= ncolLX ; i++){
    LambdaX->p[i] = pLXP[i];
  }
                                   
  if(path){
    exLambdaX = cs_kroneckerSI(LambdaX, nlGR[nG]);
  }

  lambda[0] = cs_spalloc(nL, 1, nL, true, false); 
  mulambda = cs_spalloc(nL, 1, nL, true, false); 

  for (i = 0 ; i < nL ; i++){
    lambda[0]->i[i] = i;
    lambda[0]->x[i] = i*0.01+0.01;
    mulambda->i[i] = i;
    mulambda->x[i] = 0.0;
  }
  lambda[0]->p[0] = 0;
  lambda[0]->p[1] = nL;
  mulambda->p[0] = 0;
  mulambda->p[1] = nL;

// read in prior covaraince matrix and mean vector for structural parameters.

  pvL = cs_spalloc(nL, nL, pow(nL, 2.0-LvpP[nL*nL]), true, false);
  pmuL = cs_spalloc(nL, 1, nL, true, false);
  Lrv = cs_spalloc(nL, 1, nL, true, false);

  cnt = 0;
  for (i = 0 ; i < nL ; i++){
   pmuL->i[i] = i;
   Lrv->i[i] = i;
   pmuL->x[i] = LmupP[i];
  }
  pmuL->p[0] = 0; 
  pmuL->p[1] = nL;
  Lrv->p[0] = 0; 
  Lrv->p[1] = nL;
  
  if(LvpP[nL*nL]>0.5){   // prior for structural parameters is diagonal
    for (i = 0 ; i < nL ; i++){
     pvL->p[i] = i; 
     pvL->i[i] = i;
     pvL->x[i] = LvpP[nL*i+i];
    }
    pvL->p[nL] = nL;
  }else{
    cnt = 0;
    for (i = 0 ; i < nL; i++){
      pvB->p[i] = i*nL; 
      for (j = 0 ; j < nL; j++){
       pvL->i[cnt] = j;
       pvL->x[cnt] = LvpP[cnt];
       cnt++;
      }
    }
    pvL->p[nL] = nL*nL;
  }
  pvLS = cs_schol(1, pvL);                    // Symbolic factorisation of pvL
  pvLL = cs_chol(pvL, pvLS);                  // cholesky factorisation of pvL^{-1} for forming N(0, pvL)
}


//Rprintf("read Z");

/*******************************************/ 
/* read in random-effects design matrix Z  */
/*******************************************/
                                                    /*     Za_1 0   0   Zb_1 0   0  */
Z = cs_spalloc(nrowZ, ncolZ, nzmaxZ, true, false);  /* Z =  0  Za_2 0    0  Zb_2 0  */
                                                    /*     0   0  Za_3  0   0  Zb_3 */
for (i = 0 ; i < nzmaxZ ; i++){
  Z->i[i] = iZP[i];
  Z->x[i] = xZP[i];
}
for (i = 0 ; i <= ncolZ ; i++){
  Z->p[i] = pZP[i];
}

/********************/ 
/* form W = [X, Z]  */
/********************/


if(nG>0){
  W = cs_cbind(X,Z);
}else{
  W = cs_spalloc(nrowX, ncolX, nzmaxX, true, false);                                                           
  for (i = 0 ; i < nzmaxX ; i++){
    W->i[i] = iXP[i];
    W->x[i] = xXP[i];
  }
  for (i = 0 ; i <= ncolX ; i++){
    W->p[i] = pXP[i];
  }
}

/*****************************************************/ 
/* read in category probability matrix for po terms  */
/*****************************************************/

if(nme>0){
  cnt4=0;
  cnt3=0;
  for (l = 0 ; l < nme ; l++){
    
    me_prior_prob[l] = cs_dense(me_prior_probP, nmeo[l], nmec[l], cnt4);
    me_post_prob[l] = cs_dense0(nmeo[l], nmec[l]);
    me_Xi[l] = cs_dense(me_XiP, nrowX, nmei[l], cnt3);

    for(i=0; i<nmeo[l]; i++){ // sample starting categories.
      
      tprob = 0.0;
      for (p = 0 ; p < nmec[l]; p++){   
        tprob += me_prior_prob[l]->x[p*nmeo[l]+i];
      }
      
      for (p = 0 ; p < nmec[l]; p++){ 
        if(tprob>0.0 && runif(0.0, 1.0)<=(me_prior_prob[l]->x[p*nmeo[l]+i]/tprob)){
          me_outcome[l*nrowX+i]=p;
          tprob = -1.0;
        }else{
          tprob -= me_prior_prob[l]->x[p*nmeo[l]+i];
        }
      }
    } // draw new category
    
    // delete design matrix enteries
    for(i=0; i<(nrowX*nmei[l]*(nmec[l]-1)); i++){
      W->x[W->p[stme[l]]+i] = 0.0;
    }
    
    // fill in design matrix
    cnt=0;
    cnt2=0;
    for(k=nG; k<nGR; k++){          // Iterate through R-structures
      dimG = GRdim[k];
      for(j=0; j<nlGR[k]; j++){     // Iterate through levels
        
        me_rterm = me_rtermP[nrterms*l+cnt];               // which me_obs does this give information on
        
        for(i=0; i<dimG; i++){      // Iterate through fixed levels
          record=cnt2+nlGR[k]*i+j;
          if(me_outcome[l*nrowX+me_rterm]!=0){
            for (p = 0 ; p < nmei[l]; p++){
           //   W->x[W->p[stme[l]+p*nmei[l]+me_outcome[l*nrowX+me_rterm]-1]+record] = me_Xi[l]->x[nrowX*p+record];
              W->x[W->p[stme[l]+(me_outcome[l*nrowX+me_rterm]-1)*nmei[l]+p]+record] = me_Xi[l]->x[nrowX*p+record];
            }
          }
        }
        cnt++;
      }
      cnt2+=nlGR[k]*dimG;
    }
    cnt3 += nrowX*nmei[l];
    cnt4 += nmeo[l]*nmec[l];
  } 
} 

Wt = cs_transpose(W, true);

/*************************************************************/
/* read in inverse numerator matrix A and pedigree/phylogeny */
/*************************************************************/

for(k=0; k<nG; k++){

  if(AtermP[k]>=0){

    dimG = GRdim[k];
    A[k] = cs_spalloc(nlGR[k], nlGR[k], nzmaxP[3+AtermP[k]], true, false);  
    bv[k] = cs_spalloc(nlGR[k], dimG, nlGR[k]*dimG, true, false);
    bv_tmp[k] = cs_spalloc(nlGR[k], dimG, nlGR[k]*dimG, true, false);

    cnt2=0;
    cnt3=0;

    for (i = 0; i < AtermP[k]; i++){
     cnt2 += nzmaxP[3+i];
     cnt3 += dimP[6+i]+1;
    }            

    for (i = 0; i < nzmaxP[3+AtermP[k]]; i++){
      A[k]->i[i] = iAP[i+cnt2];
      A[k]->x[i] = xAP[i+cnt2];
    }
    for (i = 0; i <= nlGR[k] ; i++){
      A[k]->p[i] = pAP[i+cnt3];
    }

// create matrix for breeding values when sampling kronecker(G,A) of dimesion ntXpedigree

    cnt=0;

    for (i = 0 ; i < dimG; i++){
      bv_tmp[k]->p[i] = cnt;
      for (j = 0 ; j < nlGR[k] ; j++){
        bv_tmp[k]->i[cnt] = j;
        bv_tmp[k]->x[cnt] = 0.0;
        cnt++;
      }
    }

    bv_tmp[k]->p[dimG] = dimG*nlGR[k];

// create matrix for breeding values when sampling G (bv'Abv) of dimesion ntXindividuals anlaysed

    cnt=0;

    for (i = 0 ; i < dimG; i++){
      bv[k]->p[i] = cnt;
      for (j = 0 ; j < nlGR[k] ; j++){
        bv[k]->i[cnt] = j;
        bv[k]->x[cnt] = 0.0;
        cnt++;
      }
    }

    bv[k]->p[dimG] = dimG*nlGR[k];
  }
}

/**********************************/ 
/* read in G/R and prior matrices */
/**********************************/
 
for (k = 0 ; k < nGR; k++){

  dimG = GRdim[k];
  Grv[k] = cs_spalloc(1, dimG, dimG, true, false);
  ivar[k] = cs_spalloc(dimG, 1, dimG, true, false);

  for (i = 0 ; i < dimG; i++){        
    Grv[k]->p[i] = i;
    Grv[k]->i[i] = 0;
  }
  Grv[k]->p[dimG] = dimG;

  if(covu>0 && k==nG){ //Random-effect residual covariances. 

    bv[k] = cs_spalloc(nlGR[k], covu, nlGR[k]*covu, true, false);

    // previous random effects in matrix form

    cnt=0;

    for (i = 0 ; i < covu; i++){
      bv[k]->p[i] = cnt;
      for (j = 0 ; j < nlGR[k] ; j++){
        bv[k]->i[cnt] = j;
        bv[k]->x[cnt] = 0.0;
        cnt++;
      }
    }

    bv[k]->p[covu] = covu*nlGR[k];

    beta_rr = cs_spalloc(dimG, covu, dimG*covu, true, false);

    cnt=0;

    for (i = 0 ; i < covu; i++){
      beta_rr->p[i] = i*dimG;     
      for (j = 0 ; j < dimG; j++){
         beta_rr->i[cnt] = j;
         beta_rr->x[cnt] = 0.0;
         cnt++;
      }
    }
    beta_rr->p[covu] = dimG*covu;

    cnt = 0;

    dimG += covu;

    muG[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
    pG[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
    Ginv_rr = cs_spalloc(dimG, dimG, dimG*dimG, true, false);

    for (i = 0 ; i < dimG; i++){
      Ginv_rr->p[i] = i*dimG;
      muG[k]->p[i] = i*dimG;
      pG[k]->p[i] = i*dimG;           
      for (j = 0 ; j < dimG; j++){
         Ginv_rr->i[cnt] = j;
         Ginv_rr->x[cnt] = GRinvP[cnt+tvc];
    		 muG[k]->i[cnt] = j;
    		 muG[k]->x[cnt] = 0.0;
         pG[k]->i[cnt] = j;
         pG[k]->x[cnt] = GRvpP[cnt+tvc];
         cnt++;
      }
    }
    muG[k]->p[dimG] = dimG*dimG;
    pG[k]->p[dimG] = dimG*dimG;
    Ginv_rr->p[dimG] = dimG*dimG;
    GinvS_rr = cs_schol(0, Ginv_rr);  
    G_rr = cs_inv(Ginv_rr);
    G[k] = cs_schur(G_rr, covu, beta_rr);
    Ginv[k] = cs_inv(G[k]);
    ldet[k] = log(cs_invR(Ginv[k], G[k]));
    ldet_rr = log(cs_invR(Ginv_rr, G_rr));
    Gtmp[k] = cs_inv(Ginv_rr);  //delete this eventually   

  }else{

    Ginv[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
    G[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
    muG[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
    pG[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);

    cnt = 0;

    for (i = 0 ; i < dimG; i++){
      Ginv[k]->p[i] = i*dimG;
      G[k]->p[i] = i*dimG;
      muG[k]->p[i] = i*dimG;
      pG[k]->p[i] = i*dimG;           
      ivar[k]->i[i] = i;
      ivar[k]->x[i] = 1.0/GRinvP[cnt+tvc+i];
      for (j = 0 ; j < dimG; j++){
        Ginv[k]->i[cnt] = j;
        G[k]->i[cnt] = j;
        muG[k]->i[cnt] = j;
        muG[k]->x[cnt] = 0.0;
        pG[k]->i[cnt] = j;
        pG[k]->x[cnt] = GRvpP[cnt+tvc];
        Ginv[k]->x[cnt] = GRinvP[cnt+tvc];
        G[k]->x[cnt] = GRinvP[cnt+tvc];
        cnt++;
      }
    }
    Ginv[k]->p[dimG] = dimG*dimG;
    G[k]->p[dimG] = dimG*dimG;
    ivar[k]->p[0] = 0;
    ivar[k]->p[1] = dimG;
    muG[k]->p[dimG] = dimG*dimG;
    pG[k]->p[dimG] = dimG*dimG;
    ldet[k] = log(cs_invR(Ginv[k], G[k]));
    Gtmp[k] = cs_inv(Ginv[k]);  //delete this eventually   
  }
  tvc += dimG*dimG;
}

tvc -= covu*covu;

cnt2=0;
cnt=0;

for (k = 0 ; k < nGR; k++){  // Antedependce priors
  dimG = nanteBP[nGR+k];
  if(nanteBP[k]!=0){
    pvAnte[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
    pmuAnte[k] = cs_spalloc(dimG, 1, dimG, true, false);
    for (i = 0 ; i < dimG; i++){
      pvAnte[k]->p[i] = i*dimG;
      pmuAnte[k]->i[i] = i;
      pmuAnte[k]->x[i] = anteBmupP[cnt];
      cnt++;
      for (j = 0 ; j < dimG; j++){
         pvAnte[k]->i[i*dimG+j] = j;
         pvAnte[k]->x[i*dimG+j] = anteBvpP[cnt2];                 
         cnt2++;
      }
    }
    pvAnte[k]->p[dimG] = dimG*dimG;
    pmuAnte[k]->p[0] = 0;
    pmuAnte[k]->p[1] = dimG;
  }
}
       
/***************************************/	
/* Read in any conditional submatrices */
/***************************************/
	
for(k = 0 ; k < nGR; k++){	
  dimG = GRdim[k];
  if(covu>0 && k==nG){
    dimG += covu;
  }
  if(((covu==0 || k!=nG) && updateP[k]!=2 && updateP[k]!=4) || (covu>0 && k==nG && updateP[nGR+1]!=2 && updateP[nGR+1]!=4)){
    CM[k] = cs_spalloc(1,1,1,true, false);
        }else{	
    cnt = 0;	
    CM[k] = cs_spalloc(dimG-splitP[k], dimG-splitP[k], (dimG-splitP[k])*(dimG-splitP[k]), true, false);
    for (i = splitP[k] ; i < dimG; i++){
      CM[k]->p[i-splitP[k]] = (i-splitP[k])*(dimG-splitP[k]);
       for (j = splitP[k] ; j < dimG; j++){
        CM[k]->i[cnt] = j-splitP[k];
        if(covu>0 && k==nG){
    		  CM[k]->x[cnt] = G_rr->x[i*dimG+j];
        }else{
    		  CM[k]->x[cnt] = G[k]->x[i*dimG+j];
        }
        cnt++;
      }
    }
    CM[k]->p[dimG-splitP[k]] = (dimG-splitP[k])*(dimG-splitP[k]);
  }
}

/**********************************/ 
/* read in alpha  prior matrices */
/**********************************/
 
if(nalpha>0){

  Alphainv = cs_spalloc(nalpha, nalpha, nzmaxAVpP[0], true, false); 
  muAlpha = cs_spalloc(nalpha, 1, nalpha, true, false); 
                                      
  for (i = 0 ; i < nzmaxAVpP[0]; i++){
    Alphainv->i[i] = iAVpP[i];
    Alphainv->x[i] = xAVpP[i];
  }
  for (i = 0 ; i < nalpha ; i++){
    Alphainv->p[i] = pAVpP[i];
    muAlpha->i[i] = i;
    muAlpha->x[i] = AmupP[i];
  }

  Alphainv->p[nalpha] = pAVpP[nalpha];
  muAlpha->p[0] = 0;
  muAlpha->p[1] = nalpha;
  AlphainvS = cs_schol(0, Alphainv);                    
  AlphainvL = cs_chol(Alphainv, AlphainvS);                

  Xalpha = cs_spalloc(nrowZ, nalpha, nrowZ*nalpha, true, false); 
                                                             
  cnt=0;
  for (i = 0; i < nalpha; i++){
    Xalpha->p[i] = cnt;
    for (j = 0; j < nrowZ; j++){
      Xalpha->i[cnt] = j;
      Xalpha->x[cnt] = 0.0;
      cnt++;
    }
  }

  Xalpha->p[nalpha] = nrowZ*nalpha;
  tXalpha = cs_transpose(Xalpha, true);

  alphaastar = cs_spalloc(nalpha, 1, nalpha, true, false);
  alphalocation = cs_spalloc(nalpha, 1, nalpha, true, false);
  alphalocation_tmp = cs_spalloc(nalpha, 1, nalpha, true, false);

  for (i = 0 ; i < nalpha; i++){
     alphaastar->i[i] = i;    
     alphalocation->i[i] = i;  
     alphalocation->x[i] = 1.0;
     alphalocation_tmp->i[i] = i;    
  }
  alphaastar->p[0] = 0;
  alphaastar->p[1] = nalpha;
  alphalocation->p[0] = 0;
  alphalocation->p[1] = nalpha;
  alphalocation_tmp->p[0] = 0;
  alphalocation_tmp->p[1] = nalpha;

  alphazstar = cs_spalloc(ny, 1, ny, true, false);
  alphapred = cs_spalloc(ny, 1, ny, true, false);

  for (i = 0 ; i < ny; i++){   
     alphazstar->i[i] = i;
     alphapred->i[i] = i;
  }
  alphazstar->p[0] = 0; 
  alphazstar->p[1] = ny;
  alphapred->p[0] = 0; 
  alphapred->p[1] = ny;
}

/*********************************/
/* Read in proposal distribution */
/*********************************/

cnt2=0;
for (k = nG ; k < nGR; k++){

  dimG = GRdim[k];

  propC[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
  muC[k] = cs_spalloc(dimG, 1, dimG, true, false);
  propC[k+nGR] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
  muC[k+nGR] = cs_spalloc(dimG, 1, dimG, true, false);

  linki[k] = cs_spalloc(dimG, 1, dimG, true, false);
  linki_tmp[k] = cs_spalloc(dimG, 1, dimG, true, false);
  linki_tmp2[k] = cs_spalloc(dimG, 1, dimG, true, false);
  predi[k] = cs_spalloc(dimG, 1, dimG, true, false);         

  cnt=0;
  for (i = 0 ; i < dimG; i++){
    propC[k]->p[i] = i*dimG;
    propC[k+nGR]->p[i] = i*dimG;
    for (j = 0 ; j < dimG; j++){
      propC[k]->i[cnt] = j;
      propC[k]->x[cnt] = propCP[cnt2];
      propC[k+nGR]->i[cnt] = j;
      propC[k+nGR]->x[cnt] = propCP[cnt2];
      cnt++;
      cnt2++;
    }
  }
  propC[k]->p[dimG] = dimG*dimG;
  propC[k+nGR]->p[dimG] = dimG*dimG;
  muC[k]->p[0] = 0;
  muC[k]->p[1] = dimG;
  muC[k+nGR]->p[0] = 0;
  muC[k+nGR]->p[1] = dimG;

  linki[k]->p[0] = 0;
  linki[k]->p[1] = dimG;
  linki_tmp[k]->p[0] = 0;
  linki_tmp[k]->p[1] = dimG;
  linki_tmp2[k]->p[0] = 0;
  linki_tmp2[k]->p[1] = dimG;
  predi[k]->p[0] = 0;
  predi[k]->p[1] = dimG;

  for (i = 0 ; i < dimG; i++){
    muC[k]->i[i] = i;
    muC[k]->x[i] = 0.0;
    muC[k+nGR]->i[i] = i;
    muC[k+nGR]->x[i] = 0.0;

    linki[k]->i[i] = i;    
    linki_tmp[k]->i[i] = i;    
    linki_tmp2[k]->i[i] = i;    
    predi[k]->i[i] = i;    
  }
  propCinv[k] = cs_inv(propC[k]);
  propCinvS[k] = cs_schol(0, propCinv[k]);   
  propCinvL[k] = cs_chol(propCinv[k], propCinvS[k]);       
  propCinv[k+nGR] = cs_inv(propC[k+nGR]);
  propCinvS[k+nGR] = cs_schol(0, propCinv[k+nGR]);   
  propCinvL[k+nGR] = cs_chol(propCinv[k+nGR], propCinvS[k+nGR]);
}

// allocate vectors for pseudo-random effects z* and [0, *a] 


//Rprintf("allocate pseudo-random effects ");

zstar = cs_spalloc(ny, 1, ny, true, false);
linky = cs_spalloc(ny, 1, ny, true, false);
mupred = cs_spalloc(ny, 1, ny, true, false);
mulinky = cs_spalloc(ny, 1, ny, true, false);

astar = cs_spalloc(dimAS, 1, dimAS, true, false);
location = cs_spalloc(dimAS, 1, dimAS, true, false);
location_tmp = cs_spalloc(dimAS, 1, dimAS, true, false);

for(i = 0 ; i < ny; i++){   
  zstar->i[i] = i;
  linky->i[i] = i;
  mupred->i[i] = i;
  mupred->x[i] = 0.0;
  mulinky->i[i] = i;
  mulinky->x[i] = 0.0;
  linky->x[i] = liabP[i];             
}

zstar->p[0] = 0; 
zstar->p[1] = ny;
linky->p[0] = 0; 
linky->p[1] = ny;
mupred->p[0] = 0; 
mupred->p[1] = ny;
mulinky->p[0] = 0; 
mulinky->p[1] = ny;

if(nL>0){    

  linky_orig = cs_spalloc(ny, 1, ny, true, false);
  for (i = 0 ; i < ny; i++){   
    linky_orig->i[i] = i;
    linky_orig->x[i] = liabP[i];                        
  }
  linky_orig->p[0] = 0; 
  linky_orig->p[1] = ny; 
}

for (i = 0 ; i < dimAS; i++){
   astar->i[i] = i;    
}
astar->p[0] = 0;
astar->p[1] = dimAS;

for (i = 0 ; i < dimAS; i++){
   location->i[i] = i;    
   location_tmp->i[i] = i;  
   location->x[i] = 0.0;    
   location_tmp->x[i] = 0.0;      
}
location->p[0] = 0;
location->p[1] = dimAS;
location_tmp->p[0] = 0;
location_tmp->p[1] = dimAS;


if(nalpha>0){

  Worig = cs_transpose(Wt, true);
  cnt = ncolX;
  cnt2 = 0;
  for (k = 0; k < nGR; k++){
    dimG = GRdim[k];
    for (i = 0; i < dimG; i++){
      if(PXtermsP[k]==1){    // parameter expanded 
        for (j = W->p[cnt]; j < W->p[cnt+nlGR[k]]; j++){ 
          W->x[j] = Worig->x[j]*alphalocation->x[cnt2];     
        }
        cnt2 ++; 
      }
      cnt += nlGR[k];
    }
  } 
}


for (k = 0 ; k < nGR ; k++){
   GinvS[k] = cs_schol(0, Ginv[k]);                    // Symbolic factorisation of G
   GinvL[k] = cs_chol(Ginv[k], GinvS[k]);              // cholesky factorisation of G^{-1} for forming N(0, G \otimes I)
}
    
// form KGinv = G^{-1} \otimes I    

for (k = 0 ; k < nG ; k++){
   if(AtermP[k]>=0){
     KGinv[k] = cs_kroneckerA(Ginv[k],A[k]);         //  form kronecker(G^{-1}, A^{-1}) structure
     KGinvS[k] = cs_schol(1, KGinv[k]);              // Symbolic factorisation of kronecker(G^{-1}, A^{-1})
   }else{
     KGinv[k] = cs_kroneckerI(Ginv[k],nlGR[k]);      //  form G^{-1} structure
   }
}


for (k = nG ; k < nGR ; k++){
  KGinv[k] = cs_kroneckerI(Ginv[k],nlGR[k]);      //  form G^{-1} structure
}

KRinv = cs_directsum(KGinv, nG, nGR);

// form WtmKRinv = W^{t}%*%KRinv  and t(t(WtmKRinv)) so its ordered correctly

WtmKRinv = cs_multiply(Wt, KRinv); 
WtmKRinvtmp = cs_transpose(WtmKRinv, TRUE); 
cs_spfree(WtmKRinv);
WtmKRinv = cs_transpose(WtmKRinvtmp, TRUE); 	
cs_spfree(WtmKRinvtmp);

// form M = WtmKRinv%*%W  

M = cs_multiply(WtmKRinv, W);

// form Omega = bdiag(pvB, KGinv) 

Omega = cs_omega(KGinv, nG, pvB);

// form MME = M + Omega; mixed model equations 

MME = cs_add(M, Omega, 1.0, 1.0); 

S = cs_schol(1, MME);                            // Symbolic factorisation - only has to be done once

// If sim/rec model then:

if(nL>0){
  
 lambdaI[0] = cs_kroneckerI(lambda[0], nrowLX);
 lambdaI[1] = cs_kroneckerI(lambda[0], nrowLX); 

 Lambda_tmp[0] = cs_multiply(LambdaX, lambdaI[0]);

 I = cs_spalloc(nrowLX, nrowLX, nrowLX, true, false); // 

 for (i = 0 ; i<nrowLX; i++){
   I->p[i] = i;                                                               
   I->i[i] = i;
   I->x[i] = 1.0;
 }
 I->p[nrowLX] =nrowLX;   

 Lambda[0] = cs_add(I, Lambda_tmp[0], 1.0, -1.0);

 if(path){
   exLambda = cs_kroneckerSI(Lambda[0],nlGR[nG]);
 }

 LambdaS = cs_sqr(1, Lambda[0], 1);                            // Symbolic LU factorisation of Lambda for determinant calculation

 LambdaLU[0] = cs_lu(Lambda[0], LambdaS, DBL_EPSILON);      

 for (i = 0 ; i<nrowLX; i++){  

   for (j = LambdaLU[0]->L->p[i] ; j<(LambdaLU[0]->L->p[i+1]); j++){
      if(LambdaLU[0]->L->i[j]==i){
        diagLambdaL[i]=j;
      }
   }
 }

 for (i = 0 ; i<nrowLX; i++){           
   for (j = LambdaLU[0]->U->p[i] ; j<(LambdaLU[0]->U->p[i+1]); j++){
      if(LambdaLU[0]->U->i[j]==i){
        diagLambdaU[i]=j;
      }
   }
 }

 detLambda[0] = 0.0;
 sign_detLambda[0] = 1;

  for (i = 0 ; i<nrowLX; i++){ 
    if(LambdaLU[0]->L->x[diagLambdaL[i]]<0.0){
      sign_detLambda[0] *= -1;
    }
    if(LambdaLU[0]->U->x[diagLambdaU[i]]<0.0){
      sign_detLambda[0] *= -1;
    }
    detLambda[0] += log(fabs(LambdaLU[0]->L->x[diagLambdaL[i]]))+log(fabs(LambdaLU[0]->U->x[diagLambdaU[i]]));
  }
  
  if(path){
     detLambda[0] *= double(nlGR[nG]);
  }
  
  cs_spfree(linky);

  if(path){  
   linky = cs_multiply(exLambda, linky_orig);   
  }else{
   linky = cs_multiply(Lambda[0], linky_orig);  
  }
  cs_sortdv(linky);
  
  ILY = cs_spalloc(ny*nL, nL, ny*nL, true, false); // form kronecker(I,y) where I has dimension nL
  
  for (i = 0 ; i<nL ; i++){
   ILY->p[i] = ny*i;                                                               
   for (j = 0 ; j < ny ; j++){
     ILY->i[j+ny*i] = j+ny*i;
     ILY->x[j+ny*i] = linky_orig->x[j];
   }
  }
  ILY->p[nL] = nL*ny;   
  
  if(path){
   Y = cs_multiply(exLambdaX, ILY);
  }else{
   Y = cs_multiply(LambdaX, ILY);
  }
  
  tY = cs_transpose(Y, true);
  
  tYKrinv = cs_multiply(tY, KRinv);
  tYKrinvY = cs_multiply(tYKrinv,Y);
  tYKrinvYS = cs_schol(1, tYKrinvY);
  tYKrinvYL = cs_chol(tYKrinvY, tYKrinvYS);
}

if(nalpha>0){

  XtmKRinv = cs_multiply(tXalpha, KRinv); 
  XtmKRinvtmp = cs_transpose(XtmKRinv, TRUE); 
  cs_spfree(XtmKRinv);
  XtmKRinv = cs_transpose(XtmKRinvtmp, TRUE); 	// double transposes so there row ordered within columns!
  cs_spfree(XtmKRinvtmp);

// form alphaM = t(alphaX)%*%solve(R)%*%alphaX%*%

  alphaM = cs_multiply(XtmKRinv, Xalpha);

// form MME = M + Omega; mixed model equations 

  alphaMME = cs_add(alphaM, Alphainv, 1.0, 1.0); 

  alphaS = cs_schol(1, alphaMME);                            // Symbolic factorisation - only has to be done once
}
	
GetRNGstate();                                   // get seed for random number generation

//Rprintf("start MCMC");

/**********************************************************************************************************/
/***** MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC *****/
/**********************************************************************************************************/

for (itt = 0; itt < (nitt+DICP[0]); itt++){

/***************************/
/* form G and R structures */
/***************************/

//Rprintf("form G and R structures");

  for (k = 0 ; k < nGR; k++){                         
    if(updateP[k]>0 || (updateP[nGR+1]>0 && covu>0 && (k==nG || k==(nG-1)))){
      if(AtermP[k]>=0){
        cs_kroneckerAupdate(Ginv[k],A[k],KGinv[k]);           //  form kronecker(G^{-1}, A^{-1}) structure
      }else{
        cs_kroneckerIupdate(Ginv[k],nlGR[k],KGinv[k]);        //  form G^{-1} structure
      }
    }  
  }

/***************/
/* free memory */
/***************/

  if(itt>0){
    for (i = 0 ; i < nGR; i++){
       if(AtermP[i]>=0){
        cs_nfree(KGinvL[i]);
      }
    }
    cs_spfree(astar_tmp);
    cs_spfree(dev);
    cs_nfree(L);

    if(nalpha>0){
      cs_nfree(alphaL);
      cs_spfree(alphaastar_tmp); 
    }

    if(nL>0){
      cs_spfree(w);
      cs_spfree(tYKrinvw);
    }
  }


  if(missing){	
    for (i = nG ; i < nGR; i++){
      cs_nfree(propCinvL[i]);
      cs_spfree(propCinv[i]);   
      cs_nfree(propCinvL[i+nGR]);
      cs_spfree(propCinv[i+nGR]);   
    }
  }
  cs_spfree(WtmKRinv);     
  cs_spfree(M);                               
  cs_spfree(MME);  

  if(nalpha>0){
    cs_spfree(tXalpha);
    cs_spfree(XtmKRinv);
    cs_spfree(alphaM);
    cs_spfree(alphaMME);
  }

  if(nL>0){
     cs_nfree(tYKrinvYL);     
     cs_spfree(tYKrinv);
     cs_spfree(tYKrinvY);
     if(path && missing){                
       cs_spfree(Y);     
       cs_spfree(tY);
     } 
  }

/********************************************/
/* update design matrices if me terms exist */
/********************************************/

if(itt>0){

  if(nme>0){
 
    for(l = 0 ; l < nme ; l++){     // Iterate over me_terms
      
      // Overwrite post_prob with prior_prob

      for(j=0; j<nmec[l]; j++){ 
        for(i=0; i<nmeo[l]; i++){ 
          me_post_prob[l]->x[j*nmeo[l]+i] = me_prior_prob[l]->x[j*nmeo[l]+i];
        }
      }
      
      // Multiply priors by likelihoods
      
      cnt=0; // counts r-terms
      cnt2=0;  // counts records
      
      // update posterior probabilities
      
      for(k=nG; k<nGR; k++){          // Iterate through R-structures
        dimG = GRdim[k];
        for(j=0; j<nlGR[k]; j++){     // Iterate through levels
          
          me_rterm = me_rtermP[nrterms*l+cnt];               // which me_obs does this give information on

          // remove current prediction
          for(i=0; i<dimG; i++){      // Iterate through fixed levels
            record=cnt2+nlGR[k]*i+j;
            linki[k]->x[i] = linky->x[record];    // get response vector          
            
            if(me_outcome[l*nrowX+me_rterm]!=0){
              for (p = 0 ; p < nmei[l]; p++){
                pred->x[record] -= me_Xi[l]->x[nrowX*p+record]*(location->x[stme[l]+(me_outcome[l*nrowX+me_rterm]-1)*nmei[l]+p]);
                W->x[W->p[stme[l]+(me_outcome[l*nrowX+me_rterm]-1)*nmei[l]+p]+record] = 0.0;
              }
            }
          }

          // get likelihoods
          
          for (cnt3 = 0 ; cnt3 < nmec[l]; cnt3++){   // Iterate over categories within me terms
            
            for(i=0; i<dimG; i++){      // Iterate through fixed levels
              record=cnt2+nlGR[k]*i+j;
              predi[k]->x[i] = pred->x[record];
              
              if(cnt3>0){
                for (p = 0 ; p < nmei[l]; p++){
                  predi[k]->x[i] += me_Xi[l]->x[nrowX*p+record]*(location->x[stme[l]+(cnt3-1)*nmei[l]+p]);
                }
              }
            }
            me_post_prob[l]->x[cnt3*nmeo[l]+me_rterm] *= exp(cs_dmvnorm(linki[k], predi[k], ldet[k], Ginv[k]));
          }
          cnt++;
        }
        cnt2+=nlGR[k]*dimG;
      }  
      
      // Sample States

      for(i=0; i<nmeo[l]; i++){ 
        tprob = 0.0;
        for (p = 0 ; p < nmec[l]; p++){   
          tprob += me_post_prob[l]->x[p*nmeo[l]+i];
        }
        for (p = 0 ; p < nmec[l]; p++){ 
          if(tprob>0.0 && runif(0.0, 1.0)<=(me_post_prob[l]->x[p*nmeo[l]+i]/tprob)){
            me_outcome[l*nrowX+i]=p;
            tprob = -1.0;
          }else{
            tprob -= me_post_prob[l]->x[p*nmeo[l]+i];
          }
        }
      } // draw new category
      
      // update W and pred
      cnt=0;
      cnt2=0;  
      for(k=nG; k<nGR; k++){          // Iterate through R-structures
        dimG = GRdim[k];
        for(j=0; j<nlGR[k]; j++){     // Iterate through levels
          
          me_rterm = me_rtermP[nrterms*l+cnt];               // which me_obs does this give information on
          
          for(i=0; i<dimG; i++){      // Iterate through fixed levels
            record=cnt2+nlGR[k]*i+j;
            if(me_outcome[l*nrowX+me_rterm]!=0){
              for (p = 0 ; p < nmei[l]; p++){
                pred->x[record] += me_Xi[l]->x[nrowX*p+record]*(location->x[stme[l]+(me_outcome[l*nrowX+me_rterm]-1)*nmei[l]+p]);
                W->x[W->p[stme[l]+(me_outcome[l*nrowX+me_rterm]-1)*nmei[l]+p]+record] = me_Xi[l]->x[nrowX*p+record];
              }
            }
          }
          cnt++;
        }
        cnt2+=nlGR[k]*dimG;
      }
    }
    cs_spfree(Wt);
    Wt = cs_transpose(W, true);
  }

  cs_spfree(pred);
  for (i = 0 ; i < nGR; i++){
    if(updateP[i]>0){
      cs_spfree(G[i]); 
    }
  }
  if(covu>0 && updateP[nGR+1]>0){
    cs_spfree(G[nG]);
    cs_spfree(G_rr);
  }
}

/******************/
/* form equations */
/******************/

//Rprintf("form equations\n");

  cs_directsumupdate(KGinv, nG, nGR, KRinv);
		
  WtmKRinv = cs_multiply(Wt, KRinv);   
		
  M = cs_multiply(WtmKRinv, W);                          // form M = WtmKRinv%*%W  
		
  cs_omegaupdate(KGinv, nG, pvB, Omega);                 // update Omega = bdiag(0, KGinv) 
		
  MME = cs_add(M, Omega, 1.0, 1.0);                      // form MME = M + Omega; mixed model equations 

/*************************/
/*  sample vectors from  */
/*      their prior      */
/*************************/
//Rprintf("sample vectors\n");	

/* beta */
		
  cnt = 0;
  
  for(i=0; i<ncolX; i++){
    Brv->x[i] = rnorm(0.0,1.0);
  }                   
  cs_ltsolve(pvBL->L, Brv->x);
  for(i=0; i<ncolX; i++){
    astar->x[i]  = Brv->x[i]+pmuB->x[i];
  }                          // sample fixed effects from their prior
  cnt = ncolX;
  
  for(k=0; k<nG; k++){

    dimG = GRdim[k];

    if(AtermP[k]>=0){        	
/* complex random effects */
      KGinvL[k] = cs_chol(KGinv[k] , KGinvS[k]);
      if(KGinvL[k]==NULL){
        PutRNGstate();
        error("G-structure %i is ill-conditioned (possibly because of ginverse): use proper priors if you haven't, or rescale data if you have\n", k+1);
      }
      for(i=0; i<(nlGR[k]*dimG); i++){
        bv_tmp[k]->x[i] = rnorm(0.0,1.0);
      }      

      cs_ltsolve(KGinvL[k]->L,  bv_tmp[k]->x);

      for(i=0; i<(nlGR[k]*dimG); i++){
        astar->x[cnt+i]  = bv_tmp[k]->x[KGinvS[k]->pinv[i]];
      }    

    }else{
// blocked random effects 			
      for(i=0; i<nlGR[k]; i++){
        for(j=0; j<dimG; j++){
          Grv[k]->x[j] = rnorm(0.0,1.0);
        }
        cs_ltsolve(GinvL[k]->L, Grv[k]->x);	
        for(j=0; j<dimG; j++){
          astar->x[j*nlGR[k]+i+cnt] = Grv[k]->x[j];
        }
      }
    }   
    cnt += dimG*nlGR[k];                                                           // form [0, a*]
  }
  
  /* residuals */
  
  cnt=0;

  for(k=nG; k<nGR; k++){

    dimG = GRdim[k];
    for(i=0; i<nlGR[k]; i++){
      for(j=0; j<dimG; j++){
        Grv[k]->x[j] = rnorm(0.0,1.0);
      }
      cs_ltsolve(GinvL[k]->L, Grv[k]->x);
      for(j=0; j<dimG; j++){
        zstar->x[j*nlGR[k]+i+cnt] = Grv[k]->x[j];
      }
    }
    cnt += dimG*nlGR[k];  
  }

/************************/
/* sample pseudo vector */
/************************/
//Rprintf("sample pseudo vectors\n");  

   cs_gaxpy(W, astar->x, zstar->x);
  
   for (i = 0 ; i < ny ; i++){
      zstar->x[i] *= -1.0;
      zstar->x[i] += linky->x[i];       // form y - z*
   }
  
   astar_tmp = cs_multiply(WtmKRinv, zstar);              // WtmKRinv(y - z*)
  
   for (i = 0 ; i < dimAS; i++){
      location->x[i] = 0.0;
   }
   for (i = 0 ; i < astar_tmp->p[astar_tmp->n]; i++){
      location->x[astar_tmp->i[i]] = astar_tmp->x[i];
   }

/*************/
/* solve MME */
/*************/
//Rprintf("solve MME\n");

   L = cs_chol(MME, S); 

//Rprintf("L formed\n");


   if(L==NULL){
     PutRNGstate();
     error("Mixed model equations singular: use a (stronger) prior\n");
   }

   for (i = 0 ; i < dimAS; i++){
      location_tmp->x[i] = 0.0;
   }

//Rprintf("start inverting\n");

   cs_ipvec (S->pinv, location->x, location_tmp->x, MME->n);	 // x = P*b 
   cs_lsolve(L->L,location_tmp->x);                                // x = L\x 
   cs_ltsolve (L->L, location_tmp->x);		                 // x = L'\x 
   cs_pvec (S->pinv, location_tmp->x, location->x, MME->n);        // b = P'*x 

//Rprintf("finish inverting\n");


   for (i = 0 ; i < dimAS ; i++){
      location->x[i] += astar->x[i];
   }

/***********************/
/* sample VCV matrices */
/***********************/
//Rprintf("sample G-VCV\n");

   pred = cs_multiply(W, location);
   cs_sortdv(pred); 
   dev = cs_add(linky, pred, 1.0, -1.0);   

   cnt2 = ncolX;

	 for(i=0; i< nG; i++){        
     dimG = GRdim[i];
     if(updateP[i]>0){            
       if(AtermP[i]>=0){        
         for(j=0; j<(dimG*nlGR[i]); j++){
           bv[i]->x[j]=location->x[cnt2+j];
         }
         tbv = cs_transpose(bv[i], true);
         bvA = cs_multiply(tbv, A[i]);
         bvAbv = cs_multiply(bvA, bv[i]);
         for(k=0; k<(dimG*dimG); k++){
           Gtmp[i]->x[k] = bvAbv->x[k] + pG[i]->x[k];
         }
         cs_spfree(tbv);               
         cs_spfree(bvA);               
         cs_spfree(bvAbv);
       }else{
         for(j=0; j<dimG; j++){
           for(k=j; k<dimG; k++){
             cnt = j*dimG+k;
             Gtmp[i]->x[cnt] = 0.0;
             for(l=0; l<nlGR[i]; l++){
               Gtmp[i]->x[cnt] += location->x[cnt2+nlGR[i]*j+l]*location->x[cnt2+nlGR[i]*k+l];
             }
             // sum-of-squares for random effect i
             Gtmp[i]->x[cnt] += pG[i]->x[cnt];
             // add prior sum-ofsquares prior$V*prior$nu (or I for corg or Diag(V) for corgh)
           }
         }
         for(j=1; j<dimG; j++){
           for(k=j; k<dimG; k++){
             Gtmp[i]->x[j+dimG*k-1] = Gtmp[i]->x[(j-1)*dimG+k];
           }
         }
       }

       switch(updateP[i]){					   
       case 1:
       	  cs_invR(Gtmp[i], Ginv[i]);
      	  G[i] = cs_rinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], GinvS[i]);	 
      	break;
      				   
      	case 2:    
      	  cs_invR(Gtmp[i], Ginv[i]);
      	  G[i] = cs_rCinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], splitP[i], CM[i]); 
      	break;	   
      				   
      	case 3:  
      	  G[i] = cs_rR(Gtmp[i], double(nlGR[i]), GRnpP[i], GinvS[i], Ginv[i], ldet[i], pG[i]);
      	break;
      
      	case 4:  
      	  G[i] = cs_rRsubinvwishart(Gtmp[i], double(nlGR[i]), splitP[i], GRnpP[i], pG[i], CM[i]);
      	break;
      
        case 5:
          G[i] = cs_rAnte(location, cnt2, dimG, nlGR[i], nanteBP[i], pmuAnte[i], pvAnte[i], A[i], AtermP[i], ivar[i]->x, nanteBP[2*nGR+i], pG[i], GRnpP[i]);
        break;
      
      	case 6:    
      	  cs_invR(Gtmp[i], Ginv[i]);
      	  G[i] = cs_rSinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], splitP[i]); 
      	break;	
     }
         
     ldet[i] = log(cs_invR(G[i], Ginv[i]));             
     cs_nfree(GinvL[i]);
     GinvL[i] = cs_chol(Ginv[i], GinvS[i]);                 // cholesky factorisation of G^{-1}. 
     
     if(GinvL[i]==NULL){
       PutRNGstate();
       error("G-structure %i is ill-conditioned: use proper priors if you haven't or rescale data if you have\n", i+1);
     }
   }
   cnt2 += nlGR[i]*dimG;  
  }

/************************/
/* Sample G-R Structure */
/************************/
//Rprintf("sample G-R-VCV\n");

  if(covu>0){

    cnt2 -= nlGR[nG-1]*GRdim[nG-1];

    if(updateP[nGR+1]>0){ 
      dimG = GRdim[nG]+covu;

      for(j=0; j<(covu*nlGR[nG]); j++){
        bv[nG]->x[j]=location->x[cnt2+j];
      }

      tbv = cs_transpose(beta_rr, true);
      bvA = cs_multiply(bv[nG],tbv);

      for(j=0; j<dimG; j++){
        for(k=j; k<dimG; k++){
          cnt = j*dimG+k;
          Gtmp[nG]->x[cnt] = 0.0;
          for(l=0; l<nlGR[nG]; l++){
            if(j<covu && k<covu){
              Gtmp[nG]->x[cnt] += location->x[cnt2+nlGR[nG]*j+l]*location->x[cnt2+nlGR[nG]*k+l];
            }
            if(j>=covu && k>=covu){
              Gtmp[nG]->x[cnt] += (dev->x[nlGR[nG]*(j-covu)+l]+bvA->x[nlGR[nG]*(j-covu)+l])*(dev->x[nlGR[nG]*(k-covu)+l]+bvA->x[nlGR[nG]*(k-covu)+l]);
            }
            if(j<covu && k>=covu){
               Gtmp[nG]->x[cnt] += location->x[cnt2+nlGR[nG]*j+l]*(dev->x[nlGR[nG]*(k-covu)+l]+bvA->x[nlGR[nG]*(k-covu)+l]);
            }
          }
          Gtmp[nG]->x[cnt] += pG[nG]->x[cnt];
        }
      }

      cs_spfree(tbv);
      cs_spfree(bvA);
      for(j=1; j<dimG; j++){
        for(k=j; k<dimG; k++){
          Gtmp[nG]->x[j+dimG*k-1] = Gtmp[nG]->x[(j-1)*dimG+k];
        }
      }
      if(diagR>0){
        cnt=0;  
        for(j=0; j<dimG; j++){
          for(k=0; k<dimG; k++){
            if(j==k){
              if(diagR==2 && cnt!=0){
                Gtmp[nG]->x[0] += Gtmp[nG]->x[cnt] - pG[nG]->x[cnt];
              }
            }else{
             Gtmp[nG]->x[cnt] = 0.0;
            }
            cnt++;
          }
        }	
      }
      
      switch(updateP[nGR+1]){

        case 1: 
          cs_invR(Gtmp[nG], Ginv_rr);
          G_rr = cs_rinvwishart(Ginv_rr, double(nlGR[nG])+GRnpP[nG], GinvS_rr);
        break;
        	
        case 2:
          cs_invR(Gtmp[nG], Ginv_rr);
          G_rr = cs_rCinvwishart(Ginv_rr, double(nlGR[nG])+GRnpP[nG], splitP[nG], CM[nG]);
        break;	   
        
        case 3:
          G_rr = cs_rR(Gtmp[nG], double(nlGR[nG]), GRnpP[nG], GinvS_rr, Ginv_rr, ldet_rr, pG[nG]);  
        break;
        
        case 4:  
          G_rr = cs_rRsubinvwishart(Gtmp[nG], double(nlGR[nG]), splitP[nG], GRnpP[nG], pG[nG], CM[nG]);
        break;
        
        case 6:    
          cs_invR(Gtmp[nG], Ginv_rr);
          G_rr = cs_rSinvwishart(Ginv_rr, double(nlGR[nG])+GRnpP[nG], splitP[nG]); 
        break;
      }
	
      if(diagR>0){
        cnt=0;  
        for(j=0; j<dimG; j++){
          for(k=0; k<dimG; k++){
      	    if(j==k){
              if(diagR==2){
                G_rr->x[cnt] = G_rr->x[0];
              }                          
            }else{
              G_rr->x[cnt] = 0.0;
            }
            cnt++;
          }
      	}	
      }

      ldet_rr = log(cs_invR(G_rr, Ginv_rr));
  
      G[nG] = cs_schur(G_rr, covu, beta_rr);

  //Form Z*
      for (j = 0; j < covu; j++){                     // iterate through kr
         for (k = 0; k < GRdim[nG]; k++){              // iterate through first ks enteries
           for(l=0; l<nlGR[nG]; l++){
             W->x[W->p[cnt2+nlGR[nG]*j+l]+k] = beta_rr->x[j*GRdim[nG]+k];
           }
         }
      }
      cs_spfree(pred);
      pred = cs_multiply(W, location);
      cs_sortdv(pred);
  
      if(updateP[nGR+1]==5){
        cs_spfree(dev);
        dev = cs_add(linky, pred, 1.0, -1.0); 
      }
  
      cs_spfree(Wt);
      Wt = cs_transpose(W, true);   
  
      ldet[nG] = log(cs_invR(G[nG], Ginv[nG])); 
      cnt=0;  
      for(j=0; j<GRdim[nG-1]; j++){
    		for(k=0; k<GRdim[nG-1]; k++){  
          G[nG-1]->x[cnt] = G_rr->x[G_rr->p[j]+k];
          cnt++; 
        }
      }
  
      ldet[nG-1] = log(cs_invR(G[nG-1], Ginv[nG-1])); 
      cs_nfree(GinvL[nG-1]);	
      GinvL[nG-1] = cs_chol(Ginv[nG-1], GinvS[nG-1]);               
      cs_nfree(GinvL[nG]);	
      GinvL[nG] = cs_chol(Ginv[nG], GinvS[nG]);    
      if(GinvL[nG]==NULL || GinvL[nG-1]==NULL){
         PutRNGstate();
         error("G-R structure is ill-conditioned: use proper priors if you haven't or rescale data if you have\n");
      }           
    }
  }

/**********************/
/* Sample R Structure */
/**********************/
//Rprintf("sample R-VCV\n");

   cnt2=0;

   for(i=nG; i<nGR; i++){  
    dimG = GRdim[i];
    if(updateP[i]>0){               
      for(j=0; j<dimG; j++){
        for(k=j; k<dimG; k++){
          cnt = j*dimG+k;
          Gtmp[i]->x[cnt] = 0.0;
          for(l=0; l<nlGR[i]; l++){
            Gtmp[i]->x[cnt] += dev->x[cnt2+nlGR[i]*j+l]*dev->x[cnt2+nlGR[i]*k+l];
          }
          Gtmp[i]->x[cnt] += pG[i]->x[cnt];
        }
      }
      for(j=1; j<dimG; j++){
        for(k=j; k<dimG; k++){
          Gtmp[i]->x[j+dimG*k-1] = Gtmp[i]->x[(j-1)*dimG+k];
        }
      }
      if(diagR>0){
        cnt=0;  
        for(j=0; j<dimG; j++){
    		  for(k=0; k<dimG; k++){
    		    if(j==k){
              if(diagR==2 && cnt!=0){
               Gtmp[i]->x[0] += Gtmp[i]->x[cnt] - pG[i]->x[cnt];
              }
            }else{
             Gtmp[i]->x[cnt] = 0.0;
            }
            cnt++;
    		  }
        }	
      }
  
      switch(updateP[i]){
  
        case 1: 
          cs_invR(Gtmp[i], Ginv[i]);
          G[i] = cs_rinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], GinvS[i]);
        break;
  				
        case 2: 
    		  cs_invR(Gtmp[i], Ginv[i]);
    		  G[i] = cs_rCinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], splitP[i], CM[i]);
    		break;	   
  	
        case 3:
         G[i] = cs_rR(Gtmp[i], double(nlGR[i]), GRnpP[i], GinvS[i], Ginv[i], ldet[i], pG[i]);  
        break;
  
    		case 4:  
    		  G[i] = cs_rRsubinvwishart(Gtmp[i], double(nlGR[i]), splitP[i], GRnpP[i], pG[i], CM[i]);
    		break;
  
        case 5:
          G[i] = cs_rAnte(dev, cnt2, dimG, nlGR[i],nanteBP[i], pmuAnte[i], pvAnte[i], A[i], AtermP[i], ivar[i]->x, nanteBP[2*nGR+i], pG[i], GRnpP[i]);
        break;
  
    		case 6:    
    		  cs_invR(Gtmp[i], Ginv[i]);
    		  G[i] = cs_rSinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], splitP[i]); 
    		break;
      }					
      if(diagR>0){
        cnt=0;  
        for(j=0; j<dimG; j++){
    		  for(k=0; k<dimG; k++){
    		    if(j==k){
              if(diagR==2){
                G[i]->x[cnt] = G[i]->x[0];
              }                          
            }else{
             G[i]->x[cnt] = 0.0;
            }
            cnt++;
          }
        }	
      }
      ldet[i] = log(cs_invR(G[i], Ginv[i]));
      cs_nfree(GinvL[i]);	
      GinvL[i] = cs_chol(Ginv[i], GinvS[i]);                 // cholesky factorisation of R^{-1} for Gibbs sampling fully missing data. 
      if(GinvL[i]==NULL){
        PutRNGstate();
        error("R-structure %i is ill-conditioned: use proper priors if you haven't or rescale data if you have\n", i+nG);
      }
    }
    cnt2 += nlGR[i]*dimG;
  }

/********************************************************/
/* update recusrive-simultaneous structural parameters  */   
/********************************************************/


  if(nL>0){     
    if(missing){
      if(path){     
        // linky_orig updated when sampling liabilities, but  Y also needs to be updated. 
        for (i = 0 ; i<nL ; i++){
          for (j = 0 ; j < ny ; j++){
            ILY->x[j+ny*i] = linky_orig->x[j];
          }
        }
        Y = cs_multiply(exLambdaX, ILY);
        tY = cs_transpose(Y, true);
      }else{      
        // sir models are tolerated if value i is missing and Lambda[i,i] = 1 and Lambda[i,\i] = 0
        // however linky_orig needs to be updated with the sampled latent variable
        for (i = 0 ; i<ny ; i++){            
          if(observedP[i]==0){
            linky_orig->x[i] = linky->x[i];
          }
        }
      }
    }
  
    w = cs_add(linky_orig, pred, 1.0, -1.0);

    tYKrinv = cs_multiply(tY, KRinv);
    tYKrinvY = cs_multiply(tYKrinv,Y);
    tYKrinvw = cs_multiply(tYKrinv,w);


// Jarrod the prior needs to be included here 
    for(i=0; i<nL; i++){
      Lrv->x[i] = rnorm(0.0,1.0);
    }


    tYKrinvYL = cs_chol(tYKrinvY, tYKrinvYS); 

    cs_ltsolve(tYKrinvYL->L, Lrv->x);       // ~ N(0, V_lambda)

    cs_ipvec (tYKrinvYS->pinv, tYKrinvw->x, w->x, nL);	 
    cs_lsolve(tYKrinvYL->L,w->x);                                
    cs_ltsolve (tYKrinvYL->L, w->x);		                 
    cs_pvec (tYKrinvYS->pinv, w->x, tYKrinvw->x, nL);       

    lambda[lambda_new] = cs_add(Lrv, lambda[lambda_old], 1.0, 1.0);

    cs_kroneckerIupdate(lambda[lambda_new], nrowLX, lambdaI[lambda_new]);  

    Lambda_tmp[lambda_new] = cs_multiply(LambdaX, lambdaI[lambda_new]);    

    Lambda[lambda_new] = cs_add(I, Lambda_tmp[lambda_new], 1.0, -1.0);

    LambdaLU[lambda_new] = cs_lu(Lambda[lambda_new], LambdaS, DBL_EPSILON);     

    detLambda[lambda_new] = 0.0;
    sign_detLambda[lambda_new] = 1;

    for (i = 0 ; i<nrowLX; i++){ 
       if(LambdaLU[lambda_new]->L->x[diagLambdaL[i]]<0.0){
         sign_detLambda[lambda_new] *= -1;
       }
       if(LambdaLU[lambda_new]->U->x[diagLambdaU[i]]<0.0){
         sign_detLambda[lambda_new] *= -1;
       }

       detLambda[lambda_new] += log(fabs(LambdaLU[lambda_new]->L->x[diagLambdaL[i]]))+log(fabs(LambdaLU[lambda_new]->U->x[diagLambdaU[i]]));
    }

    if(path){
      detLambda[lambda_new] *= double(nlGR[nG]);
    }

    log_alphaL = detLambda[lambda_new] - detLambda[lambda_old];

    lambda_dev = cs_add(lambda[lambda_new], tYKrinvw, 1.0, -1.0);

    tl = cs_transpose(lambda_dev, true);
    tlV = cs_multiply(tl, tYKrinvY);
    tlVl = cs_multiply(tlV, lambda_dev);

    log_alphaL -= 0.5*tlVl->x[0];

    cs_spfree(lambda_dev);
    cs_spfree(tl);
    cs_spfree(tlV);
    cs_spfree(tlVl);

    lambda_dev = cs_add(lambda[lambda_old], tYKrinvw, 1.0, -1.0);

    tl = cs_transpose(lambda_dev, true);
    tlV = cs_multiply(tl, tYKrinvY);
    tlVl = cs_multiply(tlV, lambda_dev);

    log_alphaL += 0.5*tlVl->x[0];

    cs_spfree(lambda_dev);
    cs_spfree(tl);
    cs_spfree(tlV);
    cs_spfree(tlVl);

    if(sign_detLambda[lambda_old]==sign_detLambda[lambda_new] && log_alphaL>log(runif(0.0,1.0))){

      cs_spfree(linky);  
      if(path){
        cs_kroneckerSIupdate(Lambda[lambda_new],nlGR[nG],exLambda);            
        linky = cs_multiply(exLambda, linky_orig);
      }else{
        linky = cs_multiply(Lambda[lambda_new], linky_orig);
      }
      cs_sortdv(linky);                               

      detLambda[0] = detLambda[lambda_new];
      detLambda[1] = detLambda[lambda_new];
      sign_detLambda[0] = sign_detLambda[lambda_new];
      sign_detLambda[1] = sign_detLambda[lambda_new];
      lambda_old =  lambda_new;
      lambda_new -= 1.0;
      lambda_new = abs(lambda_new);

    }
    cs_spfree(Lambda_tmp[lambda_new]);                                     
    cs_spfree(Lambda[lambda_new]); 
    cs_nfree(LambdaLU[lambda_new]);     
    cs_spfree(lambda[lambda_new]);  
  }

/**********************/
/* calculate deviance */   
/**********************/
//Rprintf("calculate deviance\n");

  dbar =0.0;

  if(itt>=burnin && DICP[0]==1){
    if(itt==nitt){
      for(k=nG; k<nGR; k++){
        dimG = GRdim[k];
        if(covu>0 && k==nG){

          dimG += covu;
          for(i=0; i<(dimG*dimG); i++){          
            G_rr->x[i] = muG[k]->x[i];
          }
          ldet_rr = log(cs_invR(G_rr, Ginv_rr));
          cs_spfree(G[nG]);
          G[nG] = cs_schur(G_rr, covu, beta_rr);
          ldet[nG] = log(cs_invR(G[nG], Ginv[nG]));

          for (i = 0; i < covu; i++){                     // iterate through kr
            for (j = 0; j < GRdim[nG]; j++){              // iterate through first ks entries
              for(l=0; l<nlGR[nG]; l++){
                 W->x[W->p[W->n+nlGR[nG]*(i-covu)+l]+j] = beta_rr->x[i*GRdim[nG]+j];
              }
            }
          }
            
          cs_spfree(Wt);
          Wt = cs_transpose(W, true);
                                
          cnt=0;
          for(i=0; i<GRdim[nG-1]; i++){
            for(j=0; j<GRdim[nG-1]; j++){  
              G[nG-1]->x[cnt] = G_rr->x[G_rr->p[i]+j];
              cnt++; 
            }
          }
          ldet[nG-1] = log(cs_invR(G[nG-1], Ginv[nG-1])); 
        }else{
          for(i=0; i<(dimG*dimG); i++){          
            G[k]->x[i] = muG[k]->x[i];
          }
        }
        ldet[k] = log(cs_invR(G[k], Ginv[k]));  
      }
      for(i=0; i<ny; i++){          
        linky->x[i] = mulinky->x[i];
        pred->x[i] = mupred->x[i];
      }
      if(nL>0){

         // Have to recalculate |L| for mean(lambda)

        lambda[lambda_new] = cs_kroneckerI(mulambda, 1); 

        cs_kroneckerIupdate(lambda[lambda_new], nrowLX, lambdaI[lambda_new]);  
        
        Lambda_tmp[lambda_new] = cs_multiply(LambdaX, lambdaI[lambda_new]);   

        Lambda[lambda_new] = cs_add(I, Lambda_tmp[lambda_new], 1.0, -1.0);

        LambdaLU[lambda_new] = cs_lu(Lambda[lambda_new], LambdaS, DBL_EPSILON);     
        
        detLambda[lambda_new] = 0.0;
        sign_detLambda[lambda_new] = 1;

        for (i = 0 ; i<nrowLX; i++){ 
          if(LambdaLU[lambda_new]->L->x[diagLambdaL[i]]<0.0){
            sign_detLambda[lambda_new] *= -1;
          }
          if(LambdaLU[lambda_new]->U->x[diagLambdaU[i]]<0.0){
            sign_detLambda[lambda_new] *= -1;
          }

          detLambda[lambda_new] += log(fabs(LambdaLU[lambda_new]->L->x[diagLambdaL[i]]))+log(fabs(LambdaLU[lambda_new]->U->x[diagLambdaU[i]]));
        }

        if(path){
          detLambda[lambda_new] *= double(nlGR[nG]);
        }

        lambda_old = lambda_new;
        lambda_new -= 1.0;
        lambda_new = abs(lambda_new);
        cs_spfree(Lambda_tmp[lambda_new]);                                     
        cs_spfree(Lambda[lambda_new]); 
        cs_nfree(LambdaLU[lambda_new]);  
        cs_spfree(lambda[lambda_new]);    
      }
    }         
    cnt2=0;
    for(k=nG; k<nGR; k++){          // Iterate through R-structures
      dimG = GRdim[k];
      for(j=0; j<nlGR[k]; j++){     // Iterate through levels
        nkeep=0;
        ncond=0;
        for(i=0; i<dimG; i++){       // Iterate through fixed levels
          record=cnt2+nlGR[k]*i+j;
          linki[k]->x[i] = linky->x[record];              
          predi[k]->x[i] =  pred->x[record];
          if(familyP[record]==1 && observedP[record]==1){
            keep[nkeep] = i;
            nkeep ++;
          }else{
            if(familyP[record]==20 && observedP[record]==1){ // don't condition on threshold liabilities
            }else{
              cond[ncond] = i;
              ncond ++;
            }
          }
        }
        if(nkeep>0){      // some gaussian observed traits
          if(ncond>0 || (nkeep+ncond)!=dimG){   // some non-gaussian or non-observed traits
            dbar += cs_dcmvnorm(linki[k], predi[k], G[k], keep, nkeep, cond, ncond);    // some gaussian observed
// to be made campatible with path/sir models because linki is actually Lambda%*%l so linki/predi/G need to be multiplied by solve(Lambda)
          }else{
            dbar += cs_dmvnorm(linki[k], predi[k], ldet[k], Ginv[k]);                    // all gaussian observed (this is OK for sir )
          }
        }
      }
      cnt2+=nlGR[k]*dimG;
    }
    if(nL>0){ // for sir/path models likelihood multiply by Jacobian 
      dbar += detLambda[lambda_old];
// for path models which are not fully Gaussian (or observed) this needs to be modified (determinant of the conditional matrix should be easy) 
    }
  }

/********************/
/* update cutpoints */
/********************/

     if(cp){  
       for(i=0; i<nordinal; i++){ 
         for(j=2; j<(ncutpointsP[i]-1); j++){ 
            newcutpoints[cumsum_ncutpoints[i]+j] = rtnorm(oldcutpoints[cumsum_ncutpoints[i]+j], sdcp[i], newcutpoints[cumsum_ncutpoints[i]+j-1], oldcutpoints[cumsum_ncutpoints[i]+j+1]);
         } 
       }
       cnt2=0;
       rterm=0;

       for(k=nG; k<nGR; k++){      // Iterate through R-structures
         for(i=0; i<dimG; i++){    // Iterate through first indiviual to find any ordinal variables
           dimG = GRdim[k];
           record=cnt2+nlGR[k]*i;
           if(familyP[record]==14 || familyP[record]==20){
             nthordinal = mfacP[rterm+i];
             if(familyP[record]==14){
               cutpointMHR[nthordinal] = dcutpoints(linky, yP, observedP, record,record+nlGR[k], oldcutpoints, newcutpoints, cumsum_ncutpoints[nthordinal], ncutpointsP[nthordinal], sdcp[nthordinal], 1.0);
             }else{
               cutpointMHR[nthordinal] = dcutpoints(pred, yP, observedP, record,record+nlGR[k], oldcutpoints, newcutpoints, cumsum_ncutpoints[nthordinal], ncutpointsP[nthordinal], sdcp[nthordinal], sqrt(G[k]->x[i*(dimG+1)]));
             }
             // to be made campatible with path/sir models because linky is actually Lambda%*%y so linky and predy need to be premultiplied by solve(Lambda)

             wncp[nthordinal] *= rACCEPT;
             zncp[nthordinal] *= rACCEPT;
             wncp[nthordinal] ++;

             if(cutpointMHR[nthordinal]>log(runif(0.0,1.0))){
               zncp[nthordinal] ++;
               accp[nthordinal] ++;
               for(j=2; j<(ncutpointsP[nthordinal]-1); j++){ 
                 oldcutpoints[cumsum_ncutpoints[nthordinal]+j] = newcutpoints[cumsum_ncutpoints[nthordinal]+j];
               }
               cutpoints_updated[nthordinal] = 1;
             }else{
               cutpoints_updated[nthordinal] = 0;
             } 
             if(itt<burnin){          
               sdcp[nthordinal] *= pow(qACCEPT, ((zncp[nthordinal]/wncp[nthordinal])-0.44));
             }
	         }	 
         }
         rterm += dimG;
         cnt2+=nlGR[k]*dimG;
       }
     }

/***********************/
/* sample liabilities  */   
/***********************/
//Rprintf("sample liabilities\n");	

   if(missing){
  
     cnt2=0;
     cnt=0;
     rterm=0;  // indexes the individual R-level so with the terms in the R-structure (1 1-dimensional, and 1 2 dimensional) {1} + {2, 3}
  
     for(k=nG; k<nGR; k++){      // Iterate through R-structures
  
       dimG = GRdim[k];
     	 propCinv[k] = cs_inv(propC[k]);
    	 propCinvL[k] = cs_chol(propCinv[k], propCinvS[k]); 
    	 propCinv[k+nGR] = cs_inv(propC[k+nGR]);
    	 propCinvL[k+nGR] = cs_chol(propCinv[k+nGR], propCinvS[k+nGR]);
  
       for(j=0; j<nlGR[k]; j++){     // Iterate through levels
  
         densityl1 = 0.0;
         densityl2 = 0.0;
  
         nthmnl = 0;          //variables for the multinomial
         mndenom1 = 1.0;
         mndenom2 = 1.0;
  
         p = k+proposal[cnt+j]*nGR;  // indexes proposal distributions
  
         for(i=0; i<dimG; i++){
            record=cnt2+nlGR[k]*i+j;
            linki_tmp[k]->x[i] = rnorm(0.0,y2P[record+ny*2]);
            linki[k]->x[i] = linky->x[record];
            predi[k]->x[i] =  pred->x[record];
         }
  
         if(mvtype[cnt+j]==1){         // can be Gibbsed  
           cs_ltsolve(GinvL[k]->L, linki_tmp[k]->x);  
           for(i=0; i<dimG; i++){
             linky->x[nlGR[k]*i+j+cnt2] = linki_tmp[k]->x[i]+predi[k]->x[i];
           }
           if(path){
  
             for(i=0; i<dimG; i++){
               linki_tmp2[k]->x[i] = 0.0;
               linki_tmp[k]->x[i] = linki_tmp[k]->x[i]+predi[k]->x[i];
             }
  
             cs_ipvec (LambdaLU[lambda_old]->pinv, linki_tmp[k]->x, linki_tmp2[k]->x, dimG);	 
             cs_lsolve(LambdaLU[lambda_old]->L, linki_tmp2[k]->x);                              
             cs_usolve (LambdaLU[lambda_old]->U, linki_tmp2[k]->x);		               
             cs_ipvec (LambdaS->q, linki_tmp2[k]->x, linki_tmp[k]->x, dimG);  
  
             for(i=0; i<dimG; i++){
               linky_orig->x[nlGR[k]*i+j+cnt2] = linki_tmp[k]->x[i];
             }
           }
         }
  	 
         if(mvtype[cnt+j]<0){         // has to be MHed
           cs_ltsolve(propCinvL[p]->L, linki_tmp[k]->x);           
  
           if(path){ // Solve Lambda%*%l using LU factorisation to recover original latent variables 
  
             for(i=0; i<dimG; i++){
               linki_tmp2[k]->x[i] = 0.0;
             }
  
             cs_ipvec (LambdaLU[lambda_old]->pinv, linki[k]->x, linki_tmp2[k]->x, dimG);	 
             cs_lsolve(LambdaLU[lambda_old]->L,linki_tmp2[k]->x);                            
             cs_usolve (LambdaLU[lambda_old]->U, linki_tmp2[k]->x);		              
             cs_ipvec (LambdaS->q, linki_tmp2[k]->x, linki[k]->x, dimG);      
  
             for(i=0; i<dimG; i++){
               linki_tmp2[k]->x[i] = 0.0;
             }
  
             cs_ipvec (LambdaLU[lambda_old]->pinv, linki_tmp[k]->x, linki_tmp2[k]->x, dimG);	 
             cs_lsolve(LambdaLU[lambda_old]->L, linki_tmp2[k]->x);                              
             cs_usolve (LambdaLU[lambda_old]->U, linki_tmp2[k]->x);		               
             cs_ipvec (LambdaS->q, linki_tmp2[k]->x, linki_tmp[k]->x, dimG);  
   
           }
         }
  
         for(i=0; i<dimG; i++){       // Iterate through fixed levels
  
           record=cnt2+nlGR[k]*i+j;
  
           if(mvtype[cnt+j]>0){break;} // has been Gibbsed, or is fully observed and Gaussian and therefore known
  
           linki_tmp[k]->x[i] += linki[k]->x[i];
  
           if(familyP[record]==1 && observedP[record]==1){;
             linki_tmp[k]->x[i] = linki[k]->x[i];
  
           }else{
  
             if(observedP[record]==1){
  
               switch(familyP[record]){
  
                 case 1:  /* Normal */
                 break;
  
                 case 2:  /* Posisson */
  
                   densityl1 += dpois(yP[record], exp(linki[k]->x[i]), true);
                   densityl2 += dpois(yP[record], exp(linki_tmp[k]->x[i]), true);
  
                 break;
  
                 case 3:  /* Nominal Multinomial Logit */
  
                   if(mvtype[cnt+j]==0){   // univraiate binary models can be slice sampled
                      if(yP[record]>0.5){
                       u = linky->x[record]-log1p(exp(linky->x[record]));  // needed for the deviance calculation
                       dbar += u;
                       u -= rexp(1.0);
                       linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]),u-log1p(-exp(u)),logitt);                          
                     }else{
                       u = -log1p(exp(linky->x[record])); // needed for the deviance calculation
                       dbar += u;
                       u -= rexp(1.0);
                       linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]), -logitt, log1p(-exp(u))-u);	
                     }			
                   }else{
  
                     if(truncP[0]){
                       if(linki_tmp[k]->x[i]<(-logitt)){linki_tmp[k]->x[i]=-logitt;}
                       if(linki_tmp[k]->x[i]>logitt){linki_tmp[k]->x[i]=logitt;}
                     }
  
                     mndenom1 += exp(linki[k]->x[i]);
                     mndenom2 += exp(linki_tmp[k]->x[i]);
                     densityl1 += yP[record]*linki[k]->x[i];
                     densityl2 += yP[record]*linki_tmp[k]->x[i];
  
                     if(mfacP[rterm+i]==nthmnl){ 
                       densityl1 -= y2P[record]*log(mndenom1);
                       densityl2 -= y2P[record]*log(mndenom2);
                       nthmnl = 0;
                       mndenom1 = 1.0;
                       mndenom2 = 1.0;
                     }else{
                       nthmnl++;
                     }
                   }
                 break;
   
                 case 4: /* Weibull */
                   densityl1 += dweibull(yP[record], 1.0, exp(-linki[k]->x[i]), true);
                   densityl2 += dweibull(yP[record], 1.0, exp(-linki_tmp[k]->x[i]), true);
                 break;
  
                 case 5: /* Exponential */
                   densityl1 += dexp(yP[record], exp(-linki[k]->x[i]), true);
                   densityl2 += dexp(yP[record], exp(-linki_tmp[k]->x[i]), true);
                   /* note in C dexp is parameterised via the mean, rather than through the rate as in R!*/
                 break;
  
                 case 6: /* Censored Gaussian */
  
                   if(linki_tmp[k]->x[i]<yP[record] || linki_tmp[k]->x[i]>y2P[record]){
                      interval = (linki_tmp[k]->x[i]-y2P[record])/(y2P[record]-yP[record]);
                      remainder = (linki_tmp[k]->x[i]-y2P[record])-double(int(interval))*(y2P[record]-yP[record]);
                     if(linki_tmp[k]->x[i]<yP[record]){
                       if(int(interval)%2==1){
                         linki_tmp[k]->x[i] = yP[record] - remainder;
                       }else{
                         linki_tmp[k]->x[i] = y2P[record] + remainder;
                      }
                    }else{
                       if(int(interval)%2==1){
                        linki_tmp[k]->x[i] = yP[record] + remainder;
                      }else{
                        linki_tmp[k]->x[i] = y2P[record] - remainder;
                     }
                   }
                 }
          		   break;
           
          		   case 7: /* Censored Poisson */
          	             densityl1 += log(ppois(y2P[record], exp(linki[k]->x[i]), true, false)-ppois(yP[record], exp(linki[k]->x[i]), true, false));
          		     densityl2 += log(ppois(y2P[record], exp(linki_tmp[k]->x[i]), true, false)-ppois(yP[record], exp(linki_tmp[k]->x[i]), true, false));
          		   break;
  				 
                 case 8: /* Censored Weibull */
                   densityl1 += log(pweibull(y2P[record], 1.0, exp(-linki[k]->x[i]), true, false)-pweibull(yP[record],1.0, exp(-linki[k]->x[i]), true, false));
                   densityl2 += log(pweibull(y2P[record], 1.0, exp(-linki_tmp[k]->x[i]), true, false)-pweibull(yP[record],1.0, exp(-linki_tmp[k]->x[i]), true, false));
                 break;                                            
  
                 case 9: /* Censored Exponential */
                   densityl1 += log(pexp(y2P[record], exp(-linki[k]->x[i]), true, false)-pexp(yP[record], exp(-linki[k]->x[i]), true, false));
                   densityl2 += log(pexp(y2P[record], exp(-linki_tmp[k]->x[i]), true, false)-pexp(yP[record], exp(-linki_tmp[k]->x[i]), true, false));
                 break;  
  				 
          		   case 10: /* Zero-inflated Gaussian */
          		   break;
          						 
          		   case 11: /* Zero-inflated Poisson */
  
                   if(mfacP[rterm+i]==0){
                     mndenom1 = dpois(yP[record], exp(linki[k]->x[i]), true);  
                     mndenom2 = dpois(yP[record], exp(linki_tmp[k]->x[i]), true);  
                   }else{
                      if(truncP[0]){
                        if(linki_tmp[k]->x[i]<(-logitt)){linki_tmp[k]->x[i]=-logitt;}
                        if(linki_tmp[k]->x[i]>logitt){linki_tmp[k]->x[i]=logitt;}
                      }
                			mndenom1 += log1p(-exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])));
                			mndenom2 += log1p(-exp(linki_tmp[k]->x[i])/(1.0+exp(linki_tmp[k]->x[i])));
                  if(yP[record]>0.5){
              			  mndenom1 = log(exp(mndenom1)+exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])));
              			  mndenom2 = log(exp(mndenom2)+exp(linki_tmp[k]->x[i])/(1.0+exp(linki_tmp[k]->x[i])));  
              			}
                    densityl1 += mndenom1;
                    densityl2 += mndenom2;
                    mndenom1 = 1.0;
                    mndenom2 = 1.0;
                  }
  
            		 break;
            						 
          		   case 12: /* Zero-inflated Weibull */ 
                 break;                                            
  				 
                 case 13: /* Zero-inflated Exponential */
          		   break;  
  
                 case 14: /* Ordered Multinomial Probit */
  
                 if(mvtype[cnt+j]==0){   // univraiate binary models can be slice sampled
                  
                   if(yP[record]>1.5){
                     u = pnorm(linky->x[record], 0.0, 1.0, true, true);
                     dbar += u;    // needed for the deviance calculation   
                     u -= rexp(1.0);
                     linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]), qnorm(u, 0.0, 1.0, true, true), probitt);
                   }else{
                     u = pnorm(linky->x[record], 0.0, 1.0, false, true);
                     dbar += u;   // needed for the deviance calculation
                     u -= rexp(1.0);
                     linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]), -probitt, qnorm(u, 0.0, 1.0, false, true));		
                   }

                 }else{
                   nthordinal = mfacP[rterm+i];

                   if(int(yP[record])==1 || int(yP[record])==(ncutpointsP[nthordinal]-1)){
                     if(int(yP[record])==1){
                       densityl1 += pnorm(linki[k]->x[i], 0.0, 1.0, false,true);
                       densityl2 += pnorm(linki_tmp[k]->x[i], 0.0, 1.0, false,true);
                     }else{
                       densityl1 += pnorm(oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]]-linki[k]->x[i], 0.0, 1.0, false,true);
                       densityl2 += pnorm(oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]]-linki_tmp[k]->x[i], 0.0, 1.0, false,true);
                     }
                   }else{
                     densityl1 += log(pnorm(oldcutpoints[int(yP[record])+cumsum_ncutpoints[nthordinal]]-linki[k]->x[i], 0.0, 1.0, true,false)-pnorm(oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]]-linki[k]->x[i], 0.0, 1.0, true,false));
                     densityl2 += log(pnorm(oldcutpoints[int(yP[record])+cumsum_ncutpoints[nthordinal]]-linki_tmp[k]->x[i], 0.0, 1.0, true,false)-pnorm(oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]]-linki_tmp[k]->x[i], 0.0, 1.0, true,false));
                   }
                 }
                 break;
  
                 case 15:  /* Hurdle Poisson */
  
                  if(mfacP[rterm+i]==0){
                   if(yP[record]>0.5){
                     mndenom1 = dpois(yP[record], exp(linki[k]->x[i]), true)-log1p(-exp(-exp(linki[k]->x[i])));
                     mndenom2 = dpois(yP[record], exp(linki_tmp[k]->x[i]), true)-log1p(-exp(-exp(linki_tmp[k]->x[i]))); 
                   } 
                  }else{
                   if(yP[record]>0.5){
                      mndenom1 = (linki[k]->x[i] - log1p(exp(linki[k]->x[i])));
                      mndenom2 = (linki_tmp[k]->x[i] - log1p(exp(linki_tmp[k]->x[i]))); 		       		 
                   }else{
                     mndenom1 += -linki[k]->x[i] - log1p(exp(-linki[k]->x[i]));
                  	 mndenom2 += -linki_tmp[k]->x[i] - log1p(exp(-linki_tmp[k]->x[i]));     
                   }
                   densityl1 += mndenom1;
                   densityl2 += mndenom2;
                   mndenom1 = 1.0;
                   mndenom2 = 1.0;
                  }
                 break;
  
                 case 16:  /* Zero-truncated Poisson */
  
                   densityl1 += dpois(yP[record], exp(linki[k]->x[i]), true)-log1p(-exp(-exp(linki[k]->x[i])));
                   densityl2 += dpois(yP[record], exp(linki_tmp[k]->x[i]), true)-log1p(-exp(-exp(linki_tmp[k]->x[i])));
  
                 break;
  
                 case 17:  /* Geometric */
  
                   densityl1 += linki[k]->x[i]-log1p(exp(linki[k]->x[i]))-yP[record]*(linki[k]->x[i]+log1p(exp(-linki[k]->x[i]))); 
                   densityl2 += linki_tmp[k]->x[i]-log1p(exp(linki_tmp[k]->x[i]))-yP[record]*(linki_tmp[k]->x[i]+log1p(exp(-linki_tmp[k]->x[i])));  
                   /* This is the geometric with support 0 */
  
                 break;
  
                 case 18:  /* Zero-altered Poisson */
  
                  if(mfacP[rterm+i]==0){
                     if(yP[record]>0.5){
                       mndenom1 = dpois(yP[record], exp(linki[k]->x[i]), true)-log1p(-exp(-exp(linki[k]->x[i])));
                       mndenom2 = dpois(yP[record], exp(linki_tmp[k]->x[i]), true)-log1p(-exp(-exp(linki_tmp[k]->x[i]))); 
                     } 
                   }else{
    	               if(yP[record]>0.5){
                       mndenom1 = log1p(-pexp(exp(linki[k]->x[i]), 1.0,true, false));
                   		 mndenom2 = log1p(-pexp(exp(linki_tmp[k]->x[i]),1.0, true, false));			       		 	       		 
          		       }else{
                       mndenom1 += pexp(exp(linki[k]->x[i]),1.0,true, true);
                 			 mndenom2 += pexp(exp(linki_tmp[k]->x[i]),1.0,true, true);     
                     }
                     densityl1 += mndenom1;
                     densityl2 += mndenom2;
                     mndenom1 = 1.0;
                     mndenom2 = 1.0;
                   }
  
                 break;
  
                 case 19:  /* Zero-inflated Binomial */
  
                   if(mfacP[rterm+i]==0){  // non-zero bit
                     mndenom1 = dbinom(yP[record], y2P[record], exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])), true);  
                     mndenom2 = dbinom(yP[record], y2P[record], exp(linki_tmp[k]->x[i])/(1.0+exp(linki_tmp[k]->x[i])), true);  
                   }else{
                      if(truncP[0]){
                        if(linki_tmp[k]->x[i]<(-logitt)){linki_tmp[k]->x[i]=-logitt;}
                        if(linki_tmp[k]->x[i]>logitt){linki_tmp[k]->x[i]=logitt;}
                      }
                			mndenom1 += log1p(-exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])));
                			mndenom2 += log1p(-exp(linki_tmp[k]->x[i])/(1.0+exp(linki_tmp[k]->x[i])));
                      if(yP[record]>0.5){
                			  mndenom1 = log(exp(mndenom1)+exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])));
                			  mndenom2 = log(exp(mndenom2)+exp(linki_tmp[k]->x[i])/(1.0+exp(linki_tmp[k]->x[i])));  
                			}
                      densityl1 += mndenom1;
                      densityl2 += mndenom2;
                      mndenom1 = 1.0;
                      mndenom2 = 1.0;
                   }
                 break;
  
                 case 20: /* Threshold */
                   if(cp){   // with cutpoints only Gibbs if MH step of cutpoints is accepted
                     nthordinal = mfacP[rterm+i];
                     if(DICP[0]==1){
                       dbar += pcmvnorm(predi[k], linki[k], G[k], i, oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]], oldcutpoints[int(yP[record])+cumsum_ncutpoints[nthordinal]]);    
                     } 
                     if(cutpoints_updated[nthordinal] == 1 || itt==0){
                         linky->x[record] = rtcmvnorm(predi[k], linki[k], G[k], i, oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]], oldcutpoints[int(yP[record])+cumsum_ncutpoints[nthordinal]]);
                     }
                     linki[k]->x[i] = linky->x[record];
                     linki_tmp[k]->x[i] = linky->x[record];
                   }else{
  
               // binary models can be directly Gibbsed
                     if(yP[record]>1.5){
                        if(DICP[0]==1){
                         dbar += pcmvnorm(predi[k], linki[k], G[k], i,  0.0, probitt);    
                       } 
                       linky->x[record] = rtcmvnorm(predi[k], linki[k], G[k], i, 0.0, probitt);
                     }else{ 
                       if(DICP[0]==1){
                         dbar += pcmvnorm(predi[k], linki[k], G[k], i,  -probitt, 0.0);    
                       } 
                       linky->x[record] = rtcmvnorm(predi[k], linki[k], G[k], i, -probitt, 0.0);	
                     }
                     linki[k]->x[i] = linky->x[record];
                     linki_tmp[k]->x[i] = linky->x[record];
                   }
                 break;
  
                 case 21: /* zitobit */
                    PutRNGstate();
                    error("sorry- zitobit not yet implemented\n");
                  break;
  
                 case 22: /* nzbinom */
                   
                   if(mvtype[cnt+j]==0){   // univraiate binary models can be slice sampled
                     if(yP[record]>0.5){
                       u = log1p(-pow(1.0-exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])), y2P[record])); // needed for the deviance calculation
                       dbar += u;
                       u -= rexp(1.0);
                       linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]), -logitt, log1p(-pow(1.0-exp(u)/(1.0+exp(u)), y2P[record])));	
                     }else{
                       u = y2P[record]*log1p(-exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])));  // needed for the deviance calculation
                       dbar += u;
                       u -= rexp(1.0);
                       linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]),y2P[record]*log1p(-exp(u)/(1.0+exp(u))),logitt); 
                    }			
                   }else{
                     
                     if(truncP[0]){
                       if(linki_tmp[k]->x[i]<(-logitt)){linki_tmp[k]->x[i]=-logitt;}
                       if(linki_tmp[k]->x[i]>logitt){linki_tmp[k]->x[i]=logitt;}
                     }
                     
                     if(yP[record]<0.5){
                       densityl1 += y2P[record]*log1p(-exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])));
                       densityl2 += y2P[record]*log1p(-exp(linki_tmp[k]->x[i])/(1.0+exp(linki_tmp[k]->x[i])));
                     }else{
                       densityl1 += log1p(-pow(1.0-exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])), y2P[record]));
                       densityl2 += log1p(-pow(1.0-exp(linki_tmp[k]->x[i])/(1.0+exp(linki_tmp[k]->x[i])), y2P[record]));
                     }
                   }
                   
                 break;
                   
                 case 23:  /* noncentral scaled t */
                   
                 densityl1 += dnt(yP[record], y2P[record+ny], linki[k]->x[i]/y2P[record], true)-log(y2P[record]);
                 densityl2 += dnt(yP[record], y2P[record+ny], linki_tmp[k]->x[i]/y2P[record], true)-log(y2P[record]);
                   
                 break;
                 
                 case 24:  /* mean-shifted scaled t */
                   
                  densityl1 += dt(yP[record]-linki[k]->x[i]/y2P[record], y2P[record+ny], true)-log(y2P[record]);
                  densityl2 += dt(yP[record]-linki_tmp[k]->x[i]/y2P[record], y2P[record+ny], true)-log(y2P[record]);
                 
                 break;

                 case 25:  /* Hurdle Binomial */
                   
                   if(mfacP[rterm+i]==0){  // non-zero bit
                     if(yP[record]>0.5){
                       mndenom1 = dbinom(yP[record], y2P[record], exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])), true);  
                       mndenom1 -= log1p(-dbinom(0.0, y2P[record], exp(linki[k]->x[i])/(1.0+exp(linki[k]->x[i])), false));
                       mndenom2 = dbinom(yP[record], y2P[record], exp(linki_tmp[k]->x[i])/(1.0+exp(linki_tmp[k]->x[i])), true); 
                       mndenom2 -= log1p(-dbinom(0.0, y2P[record], exp(linki_tmp[k]->x[i])/(1.0+exp(linki_tmp[k]->x[i])), false)); 
                     } 
                   }else{
                      if(truncP[0]){
                        if(linki_tmp[k]->x[i]<(-logitt)){linki_tmp[k]->x[i]=-logitt;}
                        if(linki_tmp[k]->x[i]>logitt){linki_tmp[k]->x[i]=logitt;}
                      }
                      if(yP[record]>0.5){
                        mndenom1 += linki[k]->x[i]-log1p(exp(linki[k]->x[i]));
                        mndenom2 += linki_tmp[k]->x[i]-log1p(exp(linki_tmp[k]->x[i]));
                      }else{
                        mndenom1 += -log1p(exp(linki[k]->x[i]));
                        mndenom2 += -log1p(exp(linki_tmp[k]->x[i]));
                      }
                      densityl1 += mndenom1;
                      densityl2 += mndenom2;
                      mndenom1 = 1.0;
                      mndenom2 = 1.0;
                   }
                 
                 break;

                 case 26:  /* ztmb */
                   
                   if(truncP[0]){
                     if(linki_tmp[k]->x[i]<(-logitt)){linki_tmp[k]->x[i]=-logitt;}
                     if(linki_tmp[k]->x[i]>logitt){linki_tmp[k]->x[i]=logitt;}
                   }

                   mndenom1 /= 1.0+exp(linki[k]->x[i]);
                   mndenom2 /= 1.0+exp(linki_tmp[k]->x[i]);

                   if(yP[record]>0.5){
                      densityl1 += (linki[k]->x[i] - log1p(exp(linki[k]->x[i])));
                      densityl2 += (linki_tmp[k]->x[i] - log1p(exp(linki_tmp[k]->x[i])));                
                   }else{
                     densityl1 += -linki[k]->x[i] - log1p(exp(-linki[k]->x[i]));
                     densityl2 += -linki_tmp[k]->x[i] - log1p(exp(-linki_tmp[k]->x[i]));     
                   }
                   if(mfacP[rterm+i]==nthmnl){ 
                     densityl1 -= log1p(-mndenom1);
                     densityl2 -= log1p(-mndenom2);
                     nthmnl = 0;
                     mndenom1 = 1.0;
                     mndenom2 = 1.0;
                   }else{
                     nthmnl++;
                   }
                   

                 break;

                 case 27:  /* Zero-truncated Multinomial Logit */
  
                   if(truncP[0]){
                     if(linki_tmp[k]->x[i]<(-logitt)){linki_tmp[k]->x[i]=-logitt;}
                     if(linki_tmp[k]->x[i]>logitt){linki_tmp[k]->x[i]=logitt;}
                   }

                   if(yP[record]>0.5){
                     mndenom1 += exp(linki[k]->x[i]);
                     mndenom2 += exp(linki_tmp[k]->x[i]);
                     densityl1 += yP[record]*linki[k]->x[i];
                     densityl2 += yP[record]*linki_tmp[k]->x[i];
                     present[i] = 1; 
                   }else{
                     present[i] = 0;
                   }

                   if(mfacP[rterm+i]==nthmnl){ 
                     if(y2P[record+ny]>0.5){
                       present[i+1] = 1;
                       densityl1 -= y2P[record]*log(mndenom1);
                       densityl2 -= y2P[record]*log(mndenom2);
                     }else{  // reference category is empty so do not include it in the normalisation
                       present[i+1] = 0;
                       densityl1 -= y2P[record]*log(mndenom1-1.0);
                       densityl2 -= y2P[record]*log(mndenom2-1.0);
                     }

                     densityl1 -= log(pkk_update(linki[k], y2P[record], present, nthmnl+2, i));
                     densityl2 -= log(pkk_update(linki_tmp[k], y2P[record], present, nthmnl+2, i));

                     nthmnl = 0;
                     mndenom1 = 1.0;
                     mndenom2 = 1.0;
                   }else{
                     nthmnl++;
                   }

                 break;
               }
             }
           }
         }
         dbar += densityl1;
  
         // Jarrod: pred and G are both still refering to Lambda%*%l so anything inolving them (e.g. slice sampling (when path(1,1) is fitted) and threshold models) need to be modified.
  
         if(mvtype[cnt+j]<0){
  
           if(path){  // convert l back to Lambda%*%l - shit coding sort it out!
  
             cs_spfree(linki_tmp2[k]);
             linki_tmp2[k] = cs_multiply(Lambda[lambda_old], linki[k]);   
             for(i=0; i<dimG; i++){
                linki[k]->x[linki_tmp2[k]->i[i]] = linki_tmp2[k]->x[i];
             }
             cs_spfree(linki_tmp2[k]);
             linki_tmp2[k] = cs_multiply(Lambda[lambda_old], linki_tmp[k]);   
             for(i=0; i<dimG; i++){
                linki_tmp[k]->x[linki_tmp2[k]->i[i]] = linki_tmp2[k]->x[i];
             }
             // don't need to worry about the Jacobian below because they cancel
           }

           densityl1 += cs_dmvnorm(linki[k], predi[k], ldet[k], Ginv[k]);
           densityl2 += cs_dmvnorm(linki_tmp[k], predi[k], ldet[k], Ginv[k]);

           zn[p] *= rACCEPT; 
           wn[p] *= rACCEPT;
           wn[p] ++;
  
  
           if((densityl2-densityl1)>log(runif(0.0,1.0))){
             for(i=0; i<dimG; i++){
               linky->x[nlGR[k]*i+j+cnt2] = linki_tmp[k]->x[i];
             }
             if(path){
  
               for(i=0; i<dimG; i++){
                 linki_tmp2[k]->x[i] = 0.0;
               }
  
               cs_ipvec (LambdaLU[lambda_old]->pinv, linki_tmp[k]->x, linki_tmp2[k]->x, dimG);	 
               cs_lsolve(LambdaLU[lambda_old]->L, linki_tmp2[k]->x);                              
               cs_usolve (LambdaLU[lambda_old]->U, linki_tmp2[k]->x);		               
               cs_ipvec (LambdaS->q, linki_tmp2[k]->x, linki_tmp[k]->x, dimG);  
  
               for(i=0; i<dimG; i++){
                 linky_orig->x[nlGR[k]*i+j+cnt2] = linki_tmp[k]->x[i];
               }
             }
             zn[p]++;
             Eaccl[k]++;  
           }
  
  /***************/
  /* Adaptive MH */
  /***************/
  
           if(itt<burnin && AMtuneP[k]==1){
             t[p] ++;
             for(i=0; i<dimG; i++){
               for(l=0; l<dimG; l++){
                propC[p]->x[i*dimG+l] *= (t[p]-1.0)/t[p];
                propC[p]->x[i*dimG+l] += muC[p]->x[i]*muC[p]->x[l];
               }
             }
            for(i=0; i<dimG; i++){
              muC[p]->x[i] *= (t[p]-1.0);
              muC[p]->x[i] += linky->x[nlGR[k]*i+j+cnt2];
              muC[p]->x[i] /= t[p];
            }				 
             for(i=0; i<dimG; i++){
               for(l=0; l<dimG; l++){
                 propC[p]->x[i*dimG+l] -= ((t[p]+1.0)/t[p])*muC[p]->x[i]*muC[p]->x[l];
                 propC[p]->x[i*dimG+l] += linky->x[nlGR[k]*i+j+cnt2]*linky->x[nlGR[k]*l+j+cnt2]/t[p];
               }
               propC[p]->x[i*(dimG+1)] += 0.001/t[p];
             }
           }
         }
       }
       if(itt<burnin && AMtuneP[k]==1){
         for(i=0; i<dimG; i++){
           for(l=0; l<dimG; l++){
             propC[k]->x[i*dimG+l] += (muC[k]->x[i]*muC[k]->x[l])/(t[k]);
             propC[k]->x[i*dimG+l] *= sd[k];
             if(t[k+nGR]>0.0){
               propC[k+nGR]->x[i*dimG+l] += (muC[k+nGR]->x[i]*muC[k+nGR]->x[l])/(t[k+nGR]);
               propC[k+nGR]->x[i*dimG+l] *= sd[k+nGR];
             }
           }
         }
    	   if(wn[k]>0.0){
                 sd[k] = pow(qACCEPT, ((zn[k]/wn[k])-alpha_star[k]));
    	   }	
    	   if(wn[k+nGR]>0.0){
                 sd[k+nGR] = pow(qACCEPT, ((zn[k+nGR]/wn[k+nGR])-alpha_star[k]));
    	   }	 
         zn[k] = 0.0; 
         wn[k] = 0.0;
         zn[k+nGR] = 0.0; 
         wn[k+nGR] = 0.0;
       }
       rterm += dimG;
       cnt2+=nlGR[k]*dimG;
       cnt+=nlGR[k];
     }
   }

/***********************/
/* update alpha for PX */
/***********************/

   if(nalpha>0){
  
     for (i = 0; i < nalpha*nrowZ; i++){
       Xalpha->x[i] = 0.0;
     }
  
     cnt = ncolX;
     cnt2 = 0;
     for (k = 0; k < nG; k++){
       dimG = GRdim[k];
       for (i = 0; i < dimG; i++){
         if(PXtermsP[k]==1){   // parameter expanded 
           for (l = 0; l < nlGR[k]; l++){ 
             for (j = Worig->p[cnt+l]; j < Worig->p[cnt+l+1]; j++){ 
               Xalpha->x[cnt2*nrowZ+Worig->i[j]] += Worig->x[j]*location->x[cnt+l];    
             }
           }
           cnt2 ++; 
         }
         cnt += nlGR[k];
       }
     }
  
     tXalpha = cs_transpose(Xalpha, true);
  
     XtmKRinv = cs_multiply(tXalpha, KRinv);   
  
     alphaM = cs_multiply(XtmKRinv, Xalpha);          
  
     alphaMME = cs_add(alphaM, Alphainv, 1.0, 1.0); 
  
     alphaL = cs_chol(alphaMME, alphaS); 
  
     if(alphaL==NULL){
  PutRNGstate();
       error("alpha equations singular: use a (stronger) prior for the alphas\n");
     }
  
     for(i=0; i<nalpha; i++){
        alphaastar->x[i] = rnorm(0.0,1.0);
     }               
  
     cs_ltsolve(AlphainvL->L, alphaastar->x);
     for(i=0; i<nalpha; i++){
        alphaastar->x[i]  += muAlpha->x[i];
     }      
                  
     cnt=0;
     for(k=nG; k<nGR; k++){
       dimG = GRdim[k];
       for(i=0; i<nlGR[k]; i++){
         for(j=0; j<dimG; j++){
           Grv[k]->x[j] = rnorm(0.0,1.0);
         }
         cs_ltsolve(GinvL[k]->L, Grv[k]->x);
         for(j=0; j<dimG; j++){
           alphazstar->x[j*nlGR[k]+i+cnt] = Grv[k]->x[j];     
         }
       }
       cnt += dimG*nlGR[k];  
     }  
  
  // form W%*%theta for non expanded parameter
  
     for (i = 0; i < ny; i++){
         alphapred->x[i] = linky->x[i];
     }
     for (l = 0; l < ncolX; l++){
        for (j = W->p[l]; j < W->p[l+1]; j++){ 
           alphapred->x[W->i[j]] -= W->x[j]*location->x[l]; 
        }
     }
  
     cnt = ncolX;
  
     for (k = 0; k < nG; k++){
       dimG = GRdim[k];
       for (i = 0; i < dimG; i++){
         if(PXtermsP[k]==0){   // non-parameter expanded 
           for (l = 0; l < nlGR[k]; l++){ 
             for (j = W->p[cnt+l]; j < W->p[cnt+l+1]; j++){ 
               alphapred->x[W->i[j]] -= W->x[j]*location->x[cnt+l];    
             }
           }
         }
         cnt += nlGR[k];
       }
     }
  
     cs_gaxpy(Xalpha, alphaastar->x, alphazstar->x);
  
     for (i = 0 ; i < ny ; i++){
       alphapred->x[i] -= alphazstar->x[i];       // form y - z*
     }
  
     alphaastar_tmp = cs_multiply(XtmKRinv, alphapred);                 // XtmKRinv(y - z*)
  
     for (i = 0 ; i < alphaastar_tmp->nzmax ; i++){
       alphalocation->x[alphaastar_tmp->i[i]] = alphaastar_tmp->x[i];
     }
  
     for (i = 0 ; i < nalpha; i++){
       alphalocation_tmp->x[i] = 0.0;
     }
  
     cs_ipvec (alphaS->pinv, alphalocation->x, alphalocation_tmp->x, alphaMME->n);	 
     cs_lsolve(alphaL->L,alphalocation_tmp->x);                                // x = L\x 
     cs_ltsolve (alphaL->L, alphalocation_tmp->x);		                 // x = L'\x 
     cs_pvec (alphaS->pinv, alphalocation_tmp->x, alphalocation->x, alphaMME->n);        // b = P'*x 
  
     for (i = 0 ; i < nalpha; i++){
       alphalocation->x[i] += alphaastar->x[i];
     }  
  
     cnt = ncolX;
     cnt2 = 0;
     for (k = 0; k < nG; k++){
       dimG = GRdim[k];
       for (i = 0; i < dimG; i++){
         if(PXtermsP[k]==1){    // parameter expanded 
           for (j = Worig->p[cnt]; j < Worig->p[cnt+nlGR[k]]; j++){ 
             W->x[j] = Worig->x[j]*alphalocation->x[cnt2];     
           }
           cnt2 ++; 
         }
         cnt += nlGR[k];
       }
     } 
     cs_spfree(Wt);
     Wt = cs_transpose(W, true);                               
   }
  
/***********************/
/*  store posterior   */
/***********************/

   if(itt%1000 == 0){
       if(verboseP[0]){
         Rprintf("\n                       MCMC iteration = %i\n",itt);
         for(i=nG; i<nGR; i++){
           if(nMH[i]>0){
                                              
             Rprintf("\n Acceptance ratio for liability set %i = %f\n", i+1-nG, Eaccl[i]/(nMH[i]*1000.0));
           }
         }
         if(cp){
           for(i=0; i<nordinal; i++){
             Rprintf("\n     Acceptance ratio for cutpoint set %i = %f\n", i+1, accp[i]/1000.0);
           }
         }
         R_FlushConsole();
         R_ProcessEvents();
      }
      for(i=nG; i<nGR; i++){
        if(nMH[i]>0){
          Eaccl[i] = 0.0;
        }
      }
      if(cp){
        for(i=0; i<nordinal; i++){
          accp[i] = 0.0;
        }
      }
   }

   if(itt>=burnin && DICP[0]==1){
     mdbar *= (itt-burnin);
     mdbar += dbar;
     mdbar /= (itt-burnin+1.0);
     for(i=0; i< ny; i++){
    	 mupred->x[i] *= (itt-burnin);
    	 mupred->x[i] += pred->x[i];
    	 mupred->x[i] /= (itt-burnin+1.0);
    	 mulinky->x[i] *= (itt-burnin);
    	 mulinky->x[i] += linky->x[i];
    	 mulinky->x[i] /= (itt-burnin+1.0);
     }
     for(i=0; i<nGR; i++){    
       dimG = GRdim[i];
       if(covu>0 && i==nG){
         dimG += covu;
    	   for(j=0; j<(dimG*dimG); j++){
    	     muG[i]->x[j] *= (itt-burnin);
         	     muG[i]->x[j] += G_rr->x[j];
    	     muG[i]->x[j] /= (itt-burnin+1.0);
         }
       }else{
    	   for(j=0; j<(dimG*dimG); j++){
    	     muG[i]->x[j] *= (itt-burnin);
     	     muG[i]->x[j] += G[i]->x[j];
	         muG[i]->x[j] /= (itt-burnin+1.0);
         }
       }				 
     }
     if(nL>0){
       for (i = 0 ; i < nL; i++){
    	   mulambda->x[i] *= (itt-burnin);
         mulambda->x[i] += lambda[lambda_old]->x[i];
         mulambda->x[i] /= (itt-burnin+1.0);
       }
     }
   }

   if((itt-burnin)%thin == 0 && itt>=burnin && itt!=nitt){
     if(DICP[0]==1){
       dbarP[post_cnt] = dbar;
     }
     if(cp){
       cnt =0;
       for(i=0; i<nordinal; i++){ 
         for(j=2; j<(ncutpointsP[i]-1); j++){ 
           CPP[cnt+post_cnt*ncutpoints_store] = oldcutpoints[cumsum_ncutpoints[i]+j];
           cnt++;
         }
       }
     }
     for (i = 0 ; i < ncolX ; i++){
       LocP[i+post_cnt*(ncolX+pr*ncolZ)] = location->x[i];
     }
     if(nL>0){
       for (i = 0 ; i < nL; i++){
         lambdaP[i+post_cnt*nL] = lambda[lambda_old]->x[i];        
       }
     }
     if(pr){
       if(nalpha>0){
         cnt=0;
         cnt2=0;
         for(k=0; k<nG; k++){
           dimG = GRdim[k];
           if(PXtermsP[k]==1){
             for(j=0; j<dimG; j++){
               for(i=0; i<nlGR[k]; i++){              
                 LocP[ncolX+post_cnt*dimAS+cnt] = location->x[ncolX+cnt]*alphalocation->x[cnt2];
                 cnt++;
               }
               cnt2++;
             }               
           }else{
             for(j=0; j<dimG; j++){
               for(i=0; i<nlGR[k]; i++){              
                 LocP[ncolX+post_cnt*dimAS+cnt] = location->x[ncolX+cnt];
                 cnt++;
               }
             }
           }
         }
       }else{
         for (i = ncolX ; i < dimAS ; i++){
           LocP[i+post_cnt*dimAS] = location->x[i];
         }
       }
     }
     if(pl==1){
       if(nL>0){
         for (i = 0 ; i < ny; i++){
           PLiabP[i+post_cnt*ny]  = linky_orig->x[i];
         }
       }else{
         for (i = 0 ; i < ny; i++){
           PLiabP[i+post_cnt*ny]  = linky->x[i];
         }
       }
     }

     cnt=0;
     cnt2=0;
     for(i=0; i<nGR; i++){
       dimG = GRdim[i];
       if(covu>0 && (i==nG || i==(nG-1))){
         if(i==nG){
           dimG += covu;
           for(j=0; j<(dimG*dimG); j++){
             VarP[cnt+post_cnt*tvc] = G_rr->x[j];
             cnt++;
           }
         }
       }else{
         if(PXtermsP[i]==0){
           for(j=0; j<(dimG*dimG); j++){
             VarP[cnt+post_cnt*tvc] = G[i]->x[j];
             cnt++;
           }
         }else{
           for(j=0; j<dimG; j++){              
             for(k=0; k<dimG; k++){
               VarP[cnt+post_cnt*tvc] = G[i]->x[j*dimG+k]*alphalocation->x[cnt2+j]*alphalocation->x[cnt2+k];
               cnt++;
             }
           }
           cnt2 += dimG;
         }
       }
     }
     post_cnt++;
   }
 }

/* read back proposal distribution */

 cnt2=0;
 for(k=nG; k<nGR; k++){     
   cnt=0;
   dimG = GRdim[k];
   for (j = 0 ; j < dimG*dimG; j++){
     propCP[cnt2] = propC[k]->x[cnt];
     cnt++;
     cnt2++;
   }
 }

 if(DICP[0]==1){
   mdbar *= (itt-burnin+1.0);
   mdbar -= dbar;
   mdbar /= (itt-burnin);
   dbarP[post_cnt] = mdbar;
   dbarP[post_cnt+1] = dbar;
 }

 PutRNGstate();

  cs_spfree(X);
  cs_spfree(Z);
  cs_spfree(W);
  cs_spfree(Wt);
  cs_spfree(KRinv);
  cs_spfree(WtmKRinv);
  cs_spfree(M);
  cs_spfree(Omega);
  cs_spfree(MME);
  cs_spfree(zstar);
  cs_spfree(astar);
  cs_spfree(astar_tmp);
  cs_spfree(location);
  cs_spfree(location_tmp);
  cs_spfree(linky);
  cs_spfree(pred);
  cs_spfree(mupred);
  cs_spfree(mulinky);
  cs_spfree(dev);                              
  cs_spfree(pvB);
  cs_spfree(pmuB);  
  cs_spfree(Brv); 

  if(nL>0){                               
    cs_spfree(LambdaX);
    cs_spfree(pvL);
    cs_spfree(pmuL);
    cs_spfree(Lrv);
    cs_spfree(Lambda_tmp[lambda_old]);  
    cs_spfree(Lambda[lambda_old]);   
    cs_spfree(lambda[lambda_old]);   
    cs_spfree(lambdaI[0]);   
    cs_spfree(lambdaI[1]);
    cs_spfree(mulambda);
    if(path){  
      cs_spfree(exLambda);  
      cs_spfree(exLambdaX);    
    }
    cs_spfree(I);
    cs_spfree(linky_orig);
    cs_spfree(ILY);
    cs_spfree(Y);
    cs_spfree(tY);
    cs_nfree(pvLL);
    cs_sfree(pvLS);
    cs_sfree(LambdaS);
    cs_nfree(LambdaLU[lambda_old]);
    cs_sfree(tYKrinvYS);
    cs_nfree(tYKrinvYL);
    cs_spfree(tYKrinv);                                     
    cs_spfree(tYKrinvY);   
    cs_spfree(tYKrinvw);   
    cs_spfree(w);   
  }                                                
  cs_nfree(L);
  cs_sfree(S);
  cs_nfree(pvBL);
  cs_sfree(pvBS);
  
  if(nme>0){
    for(i=0; i<nme; i++){
      cs_spfree(me_prior_prob[i]);
      cs_spfree(me_post_prob[i]);
      if(nmei[i]>0){
        cs_spfree(me_Xi[i]);
      }
    }
  }  

  for(i=0; i<nGR; i++){
    cs_spfree(Ginv[i]);
    cs_spfree(muG[i]);
    cs_spfree(G[i]);
    cs_spfree(Gtmp[i]);
    cs_spfree(Grv[i]);
    cs_spfree(pG[i]);
    cs_spfree(CM[i]);
    cs_sfree(GinvS[i]);
    cs_nfree(GinvL[i]);
    cs_spfree(KGinv[i]);
    if(AtermP[i]>=0){  
      cs_sfree(KGinvS[i]);
      cs_nfree(KGinvL[i]);
      cs_spfree(A[i]);               
      cs_spfree(bv[i]);                                                              
      cs_spfree(bv_tmp[i]);  
    }
    if(nanteBP[i]!=0){
      cs_spfree(pvAnte[i]);
      cs_spfree(pmuAnte[i]);
    }
    cs_spfree(ivar[i]);
  }

  if(nalpha>0){
    cs_spfree(Alphainv);
    cs_spfree(muAlpha);
    cs_spfree(Xalpha);
    cs_spfree(tXalpha);
    cs_spfree(XtmKRinv);
    cs_spfree(alphaM);
    cs_spfree(alphaMME);
    cs_spfree(Worig);
    cs_sfree(alphaS);
    cs_nfree(alphaL);
    cs_sfree(AlphainvS);                  
    cs_nfree(AlphainvL);    
    cs_spfree(alphaastar);
    cs_spfree(alphazstar);
    cs_spfree(alphapred);
    cs_spfree(alphaastar_tmp);
    cs_spfree(alphalocation);
    cs_spfree(alphalocation_tmp);
  }

  if(covu>0){
    cs_sfree(GinvS_rr); 
    cs_spfree(Ginv_rr);
    cs_spfree(G_rr);
    cs_spfree(beta_rr);
    cs_spfree(bv[nG]);   
  }

  for(i=nG; i<nGR; i++){
    cs_spfree(propC[i]);
    cs_spfree(propCinv[i]);
    cs_sfree(propCinvS[i]);
    cs_nfree(propCinvL[i]);
    cs_spfree(propC[i+nGR]);
    cs_spfree(propCinv[i+nGR]);
    cs_sfree(propCinvS[i+nGR]);
    cs_nfree(propCinvL[i+nGR]);
    cs_spfree(muC[i]);
    cs_spfree(muC[i+nGR]);
    cs_spfree(linki[i]);
    cs_spfree(linki_tmp[i]);
    cs_spfree(linki_tmp2[i]);
    cs_spfree(predi[i]);
  }

  delete [] cond;
  delete [] keep;
  delete [] nMH;
  delete [] Eaccl;
  delete [] alpha_star;
  delete [] t;
  delete [] sd;
  delete [] wn; 
  delete [] zn;
  delete [] ldet;
  delete [] cumsum_ncutpoints;
  delete [] oldcutpoints;
  delete [] newcutpoints;
  delete [] cutpoints_updated;
  delete [] sdcp;
  delete [] wncp;
  delete [] zncp;
  delete [] accp;
  delete [] cutpointMHR;
  delete [] Ginv;
  delete [] muG;	
  delete [] propC;
  delete [] propCinv;
  delete [] muC;
  delete [] linki;
  delete [] linki_tmp;
  delete [] linki_tmp2;
  delete [] predi;
  delete [] G;
  delete [] pG;
  delete [] CM;
  delete [] Gtmp;
  delete [] Grv; 
  delete [] GinvS;
  delete [] GinvL;
  delete [] propCinvS;
  delete [] propCinvL;
  delete [] KGinv;
  delete [] KGinvL;  
  delete [] KGinvS;
  delete [] A;
  delete [] bv;
  delete [] bv_tmp;
  delete [] diagLambdaL;
  delete [] diagLambdaU;
  delete [] lambda;
  delete [] lambdaI;
  delete [] Lambda;
  delete [] Lambda_tmp;
  delete [] LambdaLU;
  delete [] pvAnte;
  delete [] pmuAnte;
  delete [] ivar;
  delete [] me_prior_prob;
  delete [] me_post_prob;
  delete [] me_Xi;
  delete [] nmec;
  delete [] nmei;
  delete [] nmeo;
  delete [] stme;
}
}

