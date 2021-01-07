#include "MCMCglmm.h"

double pkk(double *prob, double size, int k){

  int i;

  double end[1];
  end[0] = 0.0;
  double *pend = end;


  double vector[k];
  for(i=0; i<k; i++){
     vector[i]=0.0;
  }
  double *p = vector;

  pkk_loop(0, p, k, 0, prob, size, pend);

  return (pend[0]);
}

void pkk_loop(int start, double *p, int k, int depth, double *prob, double size, double *pend){

  int i;
  for(i=start; i<k; i++){
      if(!(start<k)){break;}
      if(depth==0){
        p[depth] = prob[i];
      }else{
        p[depth] = p[depth-1]+prob[i];  
      }
      pend[0] += pow(-1.0, k-depth+1)*pow(p[depth], size);
      pkk_loop(i+1, p, k, depth+1, prob, size, pend);
  }
}

double pkk_update(const cs *linki, double size, int *present, int K, int final_i){

  int i, k, cnt;
  double sump = 0.0;

  double end[1];
  end[0] = 0.0;
  double *pend = end;

  int start_i = final_i-K+2;

  k = 0;

  for(i=start_i; i<=(final_i+1); i++){
     if(present[i]==1){k++;}
  }

  if(k==1){
    return(1.0);
  }else{

    double vector[k];
    double prob_vector[k];

    for(i=0; i<k; i++){
       vector[i]=0.0;
    }

    cnt=0;

    for(i=start_i; i<=final_i; i++){
       if(present[i]==1){  
         prob_vector[cnt] = exp(linki->x[i]);
         sump += prob_vector[cnt];
         cnt++;
       }   
    }
    if(present[final_i+1]==1){ 
        prob_vector[cnt] = 1.0;
        sump += 1.0;
    }
    for(i=0; i<k; i++){
       prob_vector[i] /= sump;
    }


    double *p = vector;
    double *prob = prob_vector;

    pkk_loop(0, p, k, 0, prob, size, pend);

    if(pend[0]<1e-16){pend[0]=1e-16;}
    
    return(pend[0]);
  }
}                
